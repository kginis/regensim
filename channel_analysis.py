import math
import csv
import time
from time import perf_counter
from pathlib import Path
import datetime

start_time = perf_counter()
# =============================================================================
# USER-EDITABLE INPUTS
# =============================================================================
identifier = 'Lynx'



# =============================================================================
# END USER-EDITABLE INPUTS
# =============================================================================

filename = Path.cwd() / f"ThermalOutput_{identifier}.csv"
#Saving it to a sheet
with filename.open("w", newline="", encoding="utf-8") as csvfile:
    print('wowo')