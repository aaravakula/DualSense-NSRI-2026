# DualSense: Dual-Biomarker Smart Wound Dressing
### NSRI Summer Research Hackathon 2026, Engineering and Technology Track

---

## Project Overview

DualSense is a conceptual engineering design for a hydrogel-based smart wound
dressing that integrates two complementary electrochemical biosensors:

1. A pH sensor using poly(2-hydroxyethyl methacrylate) hydrogel with bromocresol
   purple indicator, targeting wound exudate pH above 7.0 as an infection signal.

2. A uric acid sensor using uricase enzyme electrodes with poly(3,4-ethylenedioxythiophene)
   conducting polymer, targeting wound exudate uric acid above 400 micromol per litre.

An AND-logic dual-biomarker fusion algorithm requires both signals to exceed threshold
simultaneously before generating an infection alert, designed to minimize false positives.

**Important:** This is a conceptual design report. No combined sensor clinical
performance data are claimed. All AND-logic and OR-logic combined performance values
are model estimates from a sensitivity analysis under stated biomarker correlation
assumptions. Prospective clinical validation is required.

---

## Repository Contents

```
dualsense_simulation.py    Main simulation and sensitivity analysis code
dualsense_results.json     JSON output from the simulation run
```

---

## Simulation Code

### Requirements

```bash
pip install numpy scipy matplotlib
```

### What the code does

The simulation has four components:

**1. Patient population simulation**
Simulates biomarker readings for 2000 patients (50% infection prevalence)
drawn from normal distributions parameterized to published clinical data.

**2. Detection logic**
Applies pH-only, uric acid-only, AND-logic, and OR-logic detection rules
to the simulated population and calculates sensitivity, specificity, PPV,
NPV, and F1 score for each configuration.

**3. AND-logic decision matrix**
Generates a 4 by 4 risk probability matrix showing estimated infection
probability for each combination of pH and uric acid categories.

**4. Cost-effectiveness calculation**
Estimates potential annual savings and clinical impact based on published
wound burden statistics from Sen et al. (2025).

### Run the simulation

```bash
python dualsense_simulation.py
```

### Example output

```
DualSense Diagnostic Performance Simulation
N patients simulated: 2000
Infection prevalence: 50%

Configuration                              Sens   Spec    PPV    NPV     F1
pH sensor only                            0.960  0.961  0.961  0.959  0.961
Uric acid sensor only                     0.994  0.995  0.995  0.994  0.995
DualSense AND-logic (proposed)            0.954  1.000  1.000  0.955  0.976
DualSense OR-logic (baseline)             1.000  0.955  0.958  1.000  0.979

KEY FINDING: AND-logic achieves highest specificity,
minimizing false positives as the research question requires.
```

Note: The simulation uses independent biomarker distributions (correlation = 0),
representing the best-case scenario for AND-logic specificity. The sensitivity
analysis in the paper models performance across correlation values from 0 to 1.

---

## Sensitivity Analysis Model

The correlated binary test model used in this paper estimates combined
AND-logic and OR-logic performance as follows:

For AND-logic sensitivity:
joint_sensitivity = (S1 x S2) + rho x sqrt(S1 x (1 minus S1) x S2 x (1 minus S2))

For AND-logic specificity:
joint_fp_rate = (FP1 x FP2) + rho x sqrt(FP1 x (1 minus FP1) x FP2 x (1 minus FP2))
joint_specificity = 1 minus joint_fp_rate

Where S1 and S2 are single-sensor sensitivities, FP1 and FP2 are single-sensor
false positive rates (1 minus specificity), and rho is the assumed biomarker
correlation coefficient ranging from 0 to 1.

The plausible biomarker correlation range is assumed as rho = 0.2 to 0.4,
based on the biological independence of bacterial alkalinization (pH) and
host immune neutrophil response (uric acid). This assumption requires
prospective validation.

---

## Key Parameters (from published literature)

| Parameter | Value | Source |
|-----------|-------|--------|
| pH infection threshold | Above 7.0 | Wang et al. (2024), n=106 DFU patients |
| pH high infection threshold | Above 8.0 | Scientific Reports (2026) multicenter |
| UA alert threshold | Above 400 micromol/L | PMC OECT bandaid (2023), bench only |
| UA high infection threshold | Above 600 micromol/L | WoundMx npj (2025), bench only |
| pH sensor sensitivity | 74 to 78% | Godau et al. (2025); ScienceDirect Janus (2025) |
| pH sensor specificity | 81 to 84% | Godau et al. (2025); ScienceDirect Janus (2025) |
| UA sensor sensitivity | 69 to 72% | PMC OECT (2023); WoundMx (2025) |
| UA sensor specificity | 85 to 88% | PMC OECT (2023); WoundMx (2025) |
| Clinical visual sensitivity | 52% | Edwards et al. (2024) meta-analysis |
| Clinical visual specificity | 46% | Edwards et al. (2024) meta-analysis |

---

## Design Limitations

1. No paired patient-level data exist for combined pH and uric acid wound monitoring.
2. Uric acid threshold evidence is from bench studies only; no clinical cohort validation.
3. Alkaline wound cleansers may cause pH false positives independent of infection.
4. Electrode drift over 24 to 72 hours must be characterized before prototyping.
5. Uricase enzyme denaturation at elevated wound temperatures is an unresolved risk.

---

## Citation

If referencing this work:
DualSense: A Conceptual Dual-Biomarker Hydrogel Wound Dressing Integrating
Electrochemical pH and Uric Acid Biosensors for Earlier Detection of Chronic
Wound Infection. NSRI Summer Research Hackathon 2026, Engineering and Technology Track.

---

## License

This project is submitted for the NSRI Summer Research Hackathon 2026.
Code is provided for transparency and reproducibility under the MIT License.
