"""Batch sweep + "most optimum" selection for regensim.

Everything is driven by the [sweep], [constraints] and [batch] tables of the
input file.  Runs are spread over worker processes (the model is pure Python,
so processes -- not threads -- are what actually use multiple cores).
"""
import csv
import itertools
import math
import multiprocessing as mp
import os
import random
import time
from pathlib import Path

import engine as E

METRIC_COLS = [
    "thrust_N", "isp_s", "Pc_bar", "min_margin", "min_rib_m", "min_width_m",
    "max_width_m", "min_channel_area_m2", "mean_channel_area_m2", "max_Twg_K",
    "max_Twl_K", "dP_bar", "coolant_dT_K", "max_coolant_velocity_m_s", "coolant_boils",
]


# ----------------------------------------------------------------------------
def expand(spec, name):
    """list | {start, stop, num|step} | scalar | "max"  ->  list of values."""
    if isinstance(spec, dict):
        a, b = spec["start"], spec["stop"]
        if "num" in spec:
            num = int(spec["num"])
            vals = [a] if num == 1 else [a + (b - a) * i / (num - 1) for i in range(num)]
        elif "step" in spec:
            step = spec["step"]
            cnt = int(math.floor((b - a) / step + 1e-9)) + 1
            vals = [a + i * step for i in range(cnt)]
            if all(float(v).is_integer() for v in (a, step)):
                vals = [int(v) for v in vals]
        else:
            raise ValueError(f"sweep '{name}': need num or step")
        return vals
    if isinstance(spec, list):
        return spec
    return [spec]


_W = {}


def _init(base, nozzle, table, constraints, batch):
    _W.update(base=base, nozzle=nozzle, table=table, con=constraints, batch=batch)


def _evaluate(cfg):
    """Returns (status, metrics dict)."""
    con, batch = _W["con"], _W["batch"]
    if "min_wall" in con and cfg["Channel_Wall"] < con["min_wall"] - 1e-9:
        return "wall", {}
    try:
        gas = E.get_gas(E.finalize(cfg), _W["nozzle"])
        two_pass = cfg["Two_Pass"]
        floor = con["min_margin"] if (batch["early_exit"] and not two_pass) else None
        out = E.simulate(cfg, gas, _W["table"], margin_floor=floor, min_rib=con["min_rib"])
    except E.ModelError:
        return "error", {}
    except Exception as e:  # keep the pool alive, report the first few at the end
        return "crash", {"_err": f"{type(e).__name__}: {e}"}
    if "rejected" in out:
        return "rib", out
    if out["aborted"] or out["min_margin"] < con["min_margin"]:
        return "margin", out
    if "max_dP_bar" in con and out["dP_bar"] > con["max_dP_bar"]:
        return "dP", out
    if "max_Twg_K" in con and out["max_Twg_K"] > con["max_Twg_K"]:
        return "Twg", out
    if con.get("no_boiling") and out["coolant_boils"]:
        return "boiling", out
    return "ok", out


def _work(task):
    """task = (gas_overrides, names, combos). One gas solve, many channel sets."""
    gas_over, names, combos = task
    rows = []
    base = _W["base"]
    for combo in combos:
        cfg = dict(base)
        cfg.update(gas_over)
        cfg.update(zip(names, combo))
        status, out = _evaluate(cfg)
        rows.append((combo, status, [out.get(k) for k in METRIC_COLS], out.get("_err")))
    return gas_over, rows


# ----------------------------------------------------------------------------
def _tasks(gas_names, gas_lists, ch_names, ch_lists, mode, samples, seed, chunk):
    if mode == "grid":
        ch_all = itertools.product(*ch_lists) if ch_lists else [()]
        for gvals in itertools.product(*gas_lists):
            gover = dict(zip(gas_names, gvals))
            it = iter(itertools.product(*ch_lists) if ch_lists else [()])
            while True:
                block = list(itertools.islice(it, chunk))
                if not block:
                    break
                yield gover, ch_names, block
    else:  # random sample, grouped by gas key so the gas solve is reused
        rng = random.Random(seed)
        groups = {}
        for _ in range(samples):
            gv = tuple(rng.choice(v) for v in gas_lists)
            cv = tuple(rng.choice(v) for v in ch_lists)
            groups.setdefault(gv, []).append(cv)
        for gv, combos in groups.items():
            gover = dict(zip(gas_names, gv))
            for i in range(0, len(combos), chunk):
                yield gover, ch_names, combos[i:i + chunk]


def run_batch(cfg, raw, workers=None, out_dir=None, max_runs=None):
    sweep_raw = raw.get("sweep", {})
    if not sweep_raw:
        raise SystemExit("No [sweep] section in input file.")
    con = {"min_margin": 1.0, "min_rib": 0.0012}
    con.update(raw.get("constraints", {}))
    batch = {"workers": 0, "performance_metric": "thrust_N", "size_metric": "min_channel_area_m2",
             "performance_tolerance": 0.0, "early_exit": True, "save_all": False,
             "mode": "grid", "samples": 100000, "seed": 1}
    batch.update(raw.get("batch", {}))
    if workers is None:
        workers = batch["workers"] or os.cpu_count() or 1

    for k in sweep_raw:
        if k not in E.DEFAULTS:
            raise SystemExit(f"Unknown sweep parameter '{k}'")
    sweep = {k: expand(v, k) for k, v in sweep_raw.items()}
    gas_names = [k for k in sweep if k in E.GAS_KEYS]
    ch_names = [k for k in sweep if k not in E.GAS_KEYS]
    gas_lists = [sweep[k] for k in gas_names]
    ch_lists = [sweep[k] for k in ch_names]
    total = math.prod(len(v) for v in sweep.values())
    if batch["mode"] == "random":
        total = min(total, batch["samples"])
    if max_runs:
        total = min(total, max_runs)
    n_gas = math.prod(len(v) for v in gas_lists) if gas_lists else 1
    print(f"Sweep: {', '.join(f'{k}[{len(v)}]' for k, v in sweep.items())}")
    print(f"Total runs: {total:,}  ({n_gas} gas-side solves)  on {workers} processes")

    base = E.finalize(dict(cfg))
    base["Min_Rib"] = con["min_rib"]
    base_raw = dict(cfg)
    base_raw["Min_Rib"] = con["min_rib"]
    p_max = max([cfg["Coolant_Inlet_Pressure_Bar"]] + list(sweep.get("Coolant_Inlet_Pressure_Bar", [])))
    nozzle = E.load_nozzle(Path(cfg["_dir"]) / cfg["nozzle_file"])
    print("Building coolant property table...")
    table = E.CoolantTable(cfg["Regen_Coolant"], p_max, cache_dir=cfg["_dir"])

    names = gas_names + ch_names
    perf_i = METRIC_COLS.index(batch["performance_metric"])
    size_i = METRIC_COLS.index(batch["size_metric"])
    counts, feasible, errors = {}, [], set()
    out_dir = Path(out_dir or Path(cfg["_dir"]) / f"batch_{cfg['identifier']}")
    out_dir.mkdir(exist_ok=True)
    header = names + METRIC_COLS + ["status"]
    all_f = (out_dir / "all_runs.csv").open("w", newline="") if batch["save_all"] else None
    feas_f = (out_dir / "feasible.csv").open("w", newline="")
    all_w = csv.writer(all_f) if all_f else None
    feas_w = csv.writer(feas_f)
    feas_w.writerow(header)
    if all_w:
        all_w.writerow(header)

    # chunk size: big enough to amortise IPC, small enough to balance load
    chunk = max(1, min(500, total // (workers * 8) or 1))
    tasks = _tasks(gas_names, gas_lists, ch_names, ch_lists, batch["mode"],
                   batch["samples"], batch["seed"], chunk)
    if max_runs:
        tasks = _limit(tasks, max_runs)

    done, t0, last = 0, time.perf_counter(), 0.0
    ctx = mp.get_context("fork" if "fork" in mp.get_all_start_methods() else "spawn")
    with ctx.Pool(workers, initializer=_init, initargs=(base, nozzle, table, con, batch)) as pool:
        for gover, rows in pool.imap_unordered(_work, tasks):
            for combo, status, metrics, err in rows:
                counts[status] = counts.get(status, 0) + 1
                if err and len(errors) < 5:
                    errors.add(err)
                line = [gover[k] for k in gas_names] + list(combo) + metrics + [status]
                if all_w:
                    all_w.writerow(line)
                if status == "ok":
                    feas_w.writerow(line)
                    feasible.append(line)
            done += len(rows)
            now = time.perf_counter()
            if now - last > 2.0:
                last = now
                rate = done / (now - t0)
                print(f"\r  {done:,}/{total:,}  {rate:,.0f} runs/s  feasible {counts.get('ok', 0):,}"
                      f"  ETA {(total - done) / max(rate, 1e-9):,.0f}s   ", end="", flush=True)
    if all_f:
        all_f.close()
    feas_f.close()
    el = time.perf_counter() - t0
    print(f"\nDone: {done:,} runs in {el:.1f}s ({done / max(el, 1e-9):,.0f} runs/s)")
    print("Outcome counts:", dict(sorted(counts.items())))
    for e in errors:
        print("  crash:", e)
    import glob, shutil, tempfile
    for d in glob.glob(os.path.join(tempfile.gettempdir(), "regensim_cea_*")):
        shutil.rmtree(d, ignore_errors=True)

    if not feasible:
        print("\nNo configuration met every constraint. Loosen the sweep ranges or constraints.")
        return None

    # ------------------------------------------------ selection
    off = len(names)
    P = lambda r: r[off + perf_i]
    S = lambda r: r[off + size_i]
    best_perf = max(P(r) for r in feasible)
    tol_floor = best_perf * (1 - batch["performance_tolerance"])
    pool_ = [r for r in feasible if P(r) >= tol_floor]
    best = max(pool_, key=lambda r: (S(r), P(r)))

    ordered = sorted(feasible, key=lambda r: (-P(r), -S(r)))
    pareto, max_size = [], -1.0
    for r in ordered:
        if S(r) > max_size:
            pareto.append(r)
            max_size = S(r)
    with (out_dir / "pareto.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(pareto)
    with (out_dir / "top_candidates.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(ordered[:200])

    # best config -> reproducible input file with "max" widths resolved
    best_cfg = dict(base_raw)
    best_cfg.update(zip(names, best[:off]))
    gas = E.get_gas(E.finalize(best_cfg), nozzle)
    best_cfg = E.resolve_widths(E.finalize(best_cfg), gas)
    best_cfg["nozzle_file"] = str(Path(cfg["_dir"]) / cfg["nozzle_file"])
    best_path = out_dir / "best_config.toml"
    _write_toml(best_path, raw, best_cfg)

    print(f"\nFeasible: {len(feasible):,}   Pareto-optimal: {len(pareto)}")
    print(f"Best performance among feasible: {best_perf:.2f}")
    print("\nMost optimum configuration "
          f"(top {batch['performance_metric']} within {100 * batch['performance_tolerance']:g}%, "
          f"then largest {batch['size_metric']}):")
    for k, v in zip(names, best[:off]):
        print(f"  {k:28s} {v}")
    for k, v in zip(METRIC_COLS, best[off:off + len(METRIC_COLS)]):
        print(f"  -> {k:26s} {v:.6g}" if isinstance(v, float) else f"  -> {k:26s} {v}")
    print(f"\nFiles in {out_dir}/: feasible.csv, pareto.csv, top_candidates.csv, best_config.toml")
    print(f"Full detail for the best config:  python main.py {best_path}")
    return best


def _limit(tasks, max_runs):
    used = 0
    for g, n, combos in tasks:
        if used >= max_runs:
            return
        combos = combos[:max_runs - used]
        used += len(combos)
        yield g, n, combos


def _fmt(v):
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, str):
        return f'"{v}"'
    return repr(v)


def _write_toml(path, raw, cfg):
    lines = ["# Generated by batch optimisation -- run with: python main.py " + path.name, ""]
    for key, val in raw.items():
        if key in E.RESERVED_SECTIONS:
            continue
        if isinstance(val, dict):
            lines.append(f"[{key}]")
            for k in val:
                lines.append(f"{k} = {_fmt(cfg[k])}")
            lines.append("")
        else:
            lines.insert(2, f"{key} = {_fmt(cfg[key])}")
    path.write_text("\n".join(lines))
