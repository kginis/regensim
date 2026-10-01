# regen cooled thermal model

<img width="1929" height="793" alt="image" src="https://github.com/user-attachments/assets/78e0617b-49d9-4263-8948-448d63e789b1" />

Headless 1d thermal model for analysis of regenerative engines, designed for Nitrous/IPA rocket engines.

### Features:
- With a .csv input for nozzle geometry and propellant inputs, calculates isentropic flow properties throughout the nozzle
- Calculates temperature in the chamber, on the chamber wall, on the coolant wall, and in the coolant channels, and coolant velocity
- Film cooling analysis, but no real-world correlation for it, so it should not be trusted without skepticism
- Two-pass cooling :D

### Todo:
- Multiple different film cooling analysis methods (SP 8124 appendix A may be worth the time, see "Film Cooling Experiment v Analytical" C. Kirchberger, G. Schlieben, and O. J. Haidn paper
- Multiple different hg prediction methods. 
