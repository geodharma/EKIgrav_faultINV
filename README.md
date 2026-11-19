# EKIgrav_faultINV

MATLAB scripts of **Regularized Ensemble Kalman Inversion (EKI)** for estimating 2D fault geometry from gravity anomaly profiles, with ensemble-based uncertainty quantification.

This repository contains the code and data used in:

> Laby, D. A., Biswas, A., Heriyanto, M. F., Sungkono, S., Bery, A. A., Ahmed, B. M., Nur Rochman, J. P. G., Purwanto, M. S., & Utama, W. (2026). *Interpreting the Gravity Signatures of Faults: Efficient 2D Fault Geometry Estimation Using Regularized Ensemble Kalman Inversion with Sensitivity and Uncertainty Assessment.* Computers & Geosciences (submitted).

---

## 1. Requirements

| Item | Details |
|---|---|
| Software | MATLAB R2019a or later (`readmatrix` is used to read the data) |
| Toolbox | Statistics and Machine Learning Toolbox (`iqr`; `corr` in the Gazelle script) |
| Optional | `rms` in `data_sintetik.m` is part of base MATLAB in recent releases (Signal Processing Toolbox in older releases) |
| Hardware | Any desktop or laptop. The synthetic test (200 members × 300 iterations) runs in under a minute |

No compilation or installation is needed.

---

## 2. Repository contents

| File | Purpose |
|---|---|
| `EKI_invfaultm1.m` | **Main script.** Runs EKI on the noisy synthetic profile of Model 1 and prints/plots the result |
| `EnKF.m` | One regularized ensemble Kalman update of the whole ensemble (Kalman gain with regularization term `lambda`, set in the main script) |
| `initialize.m` | Draws the initial ensemble uniformly between the bounds `Xmin` and `Xmax` |
| `bound.m` | Keeps updated parameters inside `[Xmin, Xmax]` by reflection at the bounds |
| `objective_function.m` | Normalized data misfit of Rao & Biswas (2021), used to rank ensemble members |
| `sesar_Eq1.m` | Forward model, Model 1: inclined fault with finite depth extent (Eq. 1 of the paper) |
| `sesar_Eq2.m` | Forward model, Model 2: inclined fault with infinite depth extent, h → ∞ (Eq. 2) |
| `sesar_Eq3.m` | Forward model, Model 3: vertical fault (Eq. 3) |
| `data_sintetik.m` | Generates the synthetic data sets (noise-free and noisy) for the three models |
| `data_fault_cleanNEW.csv` | Noise-free synthetic data |
| `data_fault_noisyNEW.csv` | Noisy synthetic data (additive Gaussian noise, σ = 5 % of the anomaly range, `rng(20)`) |
| `Gravity fault(field_data).xlsx` | Field data: sheet `Gazelle Fault` (Gazelle Fault, Egypt) and sheet `Mesra` (Mersa Matruh, Egypt) |
| `EKI_invfault_Mersa.m` | Field case 1: Mersa Matruh fault, NW Egypt (vertical fault, `sesar_Eq3.m`) |
| `EKI_invfault_Gazelle.m` | Field case 2: Gazelle fault, Aswan, Egypt (inclined fault, `sesar_Eq1.m`) |
| `add_regional.m` | Adds a regional field (constant or linear) to the calculated anomaly, solved by least squares for each member. Used by the field scripts |

("sesar" means "fault" and "sintetik" means "synthetic" in Indonesian.)

### Model parameters

| Model | Forward function | Parameter vector `m` |
|---|---|---|
| 1. Inclined fault, finite depth extent | `sesar_Eq1` | `[K, z, h, alpha, x0]` |
| 2. Inclined fault, h → ∞ | `sesar_Eq2` | `[K, z, alpha, x0]` |
| 3. Vertical fault | `sesar_Eq3` | `[K, z, x0]` |

`K` is the amplitude coefficient (mGal); `z` and `h` are the depth parameters defined in Eqs. (1)–(3) of the paper (same length unit as `x`); `alpha` is the dip angle (degrees); `x0` is the horizontal position of the fault (same length unit as `x`).

### Data format

**Synthetic CSV files** (no header, 61 rows):

| Column 1 | Column 2 | Column 3 | Column 4 |
|---|---|---|---|
| x (m), 0–300 m every 5 m | Δg Model 1 (mGal) | Δg Model 2 (mGal) | Δg Model 3 (mGal) |

True models used to generate the data:

| Model | True parameters |
|---|---|
| 1 | K = 150, z = 15, h = 9, alpha = 25°, x0 = 180 |
| 2 | K = 100, z = 15, alpha = 45°, x0 = 180 |
| 3 | K = 100, z = 16, x0 = 90 |

**Field data (xlsx)**: each sheet has two columns: position along the profile and observed gravity anomaly (mGal). Rows containing NaN are removed by the scripts.

| Sheet | Location | Points | x unit |
|---|---|---|---|
| `Mesra` | Mersa Matruh, NW Egypt | 22 | km |
| `Gazelle Fault` | Gazelle (Gazal) fault, Aswan, Egypt | 17 | m |

---

## 3. Quick start: synthetic test (Model 1)

1. Download or clone this repository and open MATLAB.
2. Set the repository folder as the MATLAB current folder.
3. Run:

```matlab
EKI_invfaultm1
```

**What happens:**

- The noisy profile of Model 1 (column 2 of `data_fault_noisyNEW.csv`) is loaded.
- An ensemble of `Ne = 200` models is drawn uniformly within ±50 % of the true parameters.
- 300 EKI iterations are performed. The Command Window prints the best misfit at every iteration:
  ```
  EKI: Iteration = 1  Error = ...
  ...
  EKI: Iteration = 300  Error = ...
  ```
- A figure is displayed comparing the observed data (red stars) with the response of the best model (blue line).
- The Command Window prints the matrix `results` (one row per parameter `K, z, h, alpha, x0`):

  | Column 1 | Column 2 | Column 3 | Column 4 |
  |---|---|---|---|
  | true value | best-fitting member | ensemble median | ensemble inter-quartile range (uncertainty) |

**Expected result:** the misfit decreases steadily and the best-fit curve follows the observed data. The best-fitting parameters are close to the true values (within a few percent for `K` and `x0`). The script fixes the random seed (`rng(1)`), so repeated runs on the same MATLAB version give identical numbers. Change or remove the seed to see the run-to-run variability of the stochastic EKI.

---

## 4. Running Models 2 and 3

Edit `EKI_invfaultm1.m` as follows.

| Line(s) in `EKI_invfaultm1.m` | Model 2 | Model 3 |
|---|---|---|
| `gobs = data(:,2);` | `gobs = data(:,3);` | `gobs = data(:,4);` |
| True parameter vector | uncomment `m_true2` and use it in place of `m_true1` | uncomment `m_true3` and use it in place of `m_true1` |
| every call to `sesar_Eq1(...)` (three calls) | `sesar_Eq2(...)` | `sesar_Eq3(...)` |

To invert the **noise-free** data, replace `data_fault_noisyNEW.csv` with `data_fault_cleanNEW.csv`.

---

## 5. Regenerating the synthetic data

`data_sintetik.m` reproduces the synthetic data sets. It uses fixed seed `rng(20)` and noise σ = 5 % of the anomaly range for each model. It prints the noise statistics and plots the profiles. To overwrite the CSV files, uncomment the `writematrix` lines at the end of the script.

---

## 6. Field data

The two field scripts have the same structure as `EKI_invfaultm1.m` and use the same functions (`initialize.m`, `EnKF.m`, `bound.m`, `objective_function.m`). Only the input data, the search space, the forward model and the control parameters are changed. One function is added, `add_regional.m`. It removes the regional field or base level by least squares inside the forward calculation, so the regional is not an EKI parameter. Both scripts read `Gravity fault(field_data).xlsx`.

```matlab
EKI_invfault_Mersa      % Mersa Matruh fault
EKI_invfault_Gazelle    % Gazelle fault
```

| Item | Mersa Matruh (`EKI_invfault_Mersa.m`) | Gazelle (`EKI_invfault_Gazelle.m`) |
|---|---|---|
| Location | NW Egypt | Aswan, Egypt |
| Forward model | Vertical fault, `sesar_Eq3.m`, `m = [K, z, x0]` | Inclined fault, `sesar_Eq1.m`, `m = [K, z, h, alpha, x0]` |
| Units | x, z, x0 in km; K in mGal | x, z, h, x0 in m; K in mGal; alpha in degrees |
| Search space | K: 0–10 mGal; z: 0–5 km; x0: 8–16 km | K: 0–50 mGal; z: 1–1000 m; h: 1–10000 m; alpha: 1–179°; x0: profile extent ± 10 % |
| Regional (`reg_mode`) | `'const'` (base level c) | `'linear'` (a + b·x) |
| `Ne`, `Niter` | 100, 100 | 200, 1000 |
| `sd_kf` | 0.1 mGal | 0.3 mGal |
| `lambda` | 1e-1 | 1e-3 |
| Seed | `rng(1)` | `rng(1)` |
| Output file | `results_Mersa_EKI.mat` | `results_Gazelle_EKI_<reg_mode>.mat` |

The Gazelle data are absolute anomalies (about 48–73 mGal) on a regional gradient, so a linear regional is used. Set `reg_mode = 'const'` to compare the two options.

The Gazelle script also:

- removes non-converged members (misfit > 10 × median) before computing the posterior statistics;
- reports the derived horizontal positions of the fault segment at depths z and h (`x_top`, `x_bot`);
- reports the Spearman correlation between alpha and x0.

Each script prints the best misfit during the iterations and the RMSE (mGal) of the best member. It plots the observed data against the best-fit curve (including the regional) and prints `results = [best, median, IQR]`, one row per parameter.

**To apply the code to another profile:** copy one of the field scripts and change four things:

1. the data file or sheet;
2. the forward function (`sesar_Eq1`, `sesar_Eq2` or `sesar_Eq3`);
3. `Xmin` and `Xmax`;
4. `reg_mode`.

---

## 7. Main control parameters of the synthetic test (in `EKI_invfaultm1.m`)

| Variable | Default | Meaning |
|---|---|---|
| `Ne` | 200 | Ensemble size |
| `Niter` | 300 | Number of EKI iterations |
| `sd_kf` | 0.02 | Standard deviation of the observation perturbation; also sets the data covariance `Cd = sd_kf^2 * I` |
| `lambda` | 1e-3 | Regularization added to the data-space covariance before computing the Kalman gain |
| `Xmin`, `Xmax` | 0.5× and 1.5× the true model | Search bounds |

### Algorithm outline

1. Draw the initial ensemble uniformly within the bounds (`initialize.m`) and compute its forward responses and misfits.
2. At each iteration, update all members with the regularized Kalman gain, `K = Cxy (Cyy + Cd + λI)^-1`, using perturbed observations and a random step length in (0, 1) (`EnKF.m`). Reflect any parameter that leaves the bounds back inside them (`bound.m`).
3. Accept an updated member only if it lowers its misfit (`objective_function.m`).
4. After the last iteration, report the ensemble median as the point estimate and the inter-quartile range as the uncertainty. Also report the lowest-misfit member.

---

## 8. License

MIT License (see `LICENSE`).

## 9. Citation

If you use this code, please cite the paper above.

## 10. Contact

Dharma Arung Laby, Department of Geophysical Engineering, Institut Teknologi Sepuluh Nopember (ITS), Surabaya, Indonesia.
