# CiPA in Python

Python conversion of [FDA/CiPA](https://github.com/FDA/CiPA), pinned to commit
`539ff80d6616895de5d943e7bfb77b346d195db4`. The original checkout is in
`reference/CiPA`; the working Python implementation is in `cipa/`.

The port includes both equation models, hERG bootstrap/CMA-ES fitting,
Hill least-squares/MCMC fitting, fixed and uncertainty action-potential
simulation, current omission, sensitivity summaries, qNet intervals, and
ordinal TdP classification with leave-one-drug-out cross-validation.
Running the package requires neither R nor a C compiler.

## Setup

Python 3.11 or later is required; Python 3.12 was used for validation.
From this directory in PowerShell:

```powershell
# The local environment has already been created and populated.
.venv\Scripts\python.exe -m cipa --help
.venv\Scripts\python.exe -m pytest -q
```

To set up a fresh environment elsewhere:

```powershell
py -3.12 -m venv .venv
.venv\Scripts\python.exe -m pip install -e ".[test]"
```

Numba compiles the Python equation functions on first use, which can take
several seconds. Subsequent runs use its cache. Bundled input data and model
defaults are resolved relative to the package, so commands work from other
directories. Results are written relative to your current directory unless
you specify an output path.

## Quick control simulation

```powershell
.venv\Scripts\python.exe -m cipa ap --beats 2 --nsave 1
```

This creates `results/ap/fixed/control/metrics.csv` and, under `dose_0/`,
`trace.csv`, `trace.png`, `beat_metrics.csv`, `states_max_dv.csv`, and
`checkpoints.csv`. A short run checks operation; the original protocol uses
1,000 beats and analyzes the last 250, which remain the defaults:

```powershell
.venv\Scripts\python.exe -m cipa ap
```

To also save five consecutive final cycles, add `--trace-beats 5`:

```powershell
.venv\Scripts\python.exe -m cipa ap -d bepridil -x 1 --trace-beats 5
```

Each dose directory then also contains `trace_cycles.csv` and
`trace_cycles.png`. With 1,000 beats, these contain beats 996–1000, while
`trace.csv` and `trace.png` still contain the single beat selected by maximum
repolarization dv. `--trace-beats` is independent of `--nsave`, which controls
how many beats contribute metrics. The cycle CSV includes all the usual state
and current columns, plus `beat`, `beat_time` (within each beat), and
`simulation_time` (since the beginning of pacing). Its `time` column starts
at zero at the first saved beat and increases continuously across the five
cycles (0–10,000 ms at CL=2,000 ms). Shared boundary times appear only once;
the new beat's initial sample is kept. qNet still resets at each beat.
You must run at least as many beats as requested by `--trace-beats`.
Existing simulations must be rerun to capture these additional trajectories.

## Drug fitting and uncertainty propagation

For example, with bepridil:

```powershell
# hERG: resample whole cells within concentration strata, then fit.
.venv\Scripts\python.exe -m cipa bootstrap -d bepridil -n 2000
.venv\Scripts\python.exe -m cipa herg-fit -d bepridil
.venv\Scripts\python.exe -m cipa herg-fit -d bepridil -i "1-2000"
.venv\Scripts\python.exe -m cipa herg-summary -d bepridil
.venv\Scripts\python.exe -m cipa herg-sensitivity -d bepridil

# Six non-hERG currents and the additional hERG Hill curve in the input data.
.venv\Scripts\python.exe -m cipa hill-fit -d bepridil
.venv\Scripts\python.exe -m cipa hill-summary -d bepridil
.venv\Scripts\python.exe -m cipa filter-li2017

# Fixed and sampled inputs, at 1 through 10 times therapeutic concentration.
.venv\Scripts\python.exe -m cipa ap -d bepridil -x "1-10"
.venv\Scripts\python.exe -m cipa ap -d bepridil -x "1-10" -i "1-2000"
.venv\Scripts\python.exe -m cipa combine results/ap/fixed
.venv\Scripts\python.exe -m cipa combine results/ap/uncertainty -n 2000
.venv\Scripts\python.exe -m cipa qnet-ci
```

Inspect `mcmc_diagnostics.json` and the channel chain CSVs before using MCMC
samples. A failed diagnostic emits a warning; samples are still saved for
inspection. Increase burn-in and sample count and check mixing as necessary.
The original full fitting and simulation settings require substantial compute.
Use separate, disjoint sample ranges in separate processes to distribute work;
each sample writes to its own directory. Do not run overlapping ranges against
the same output directory.

`filter-li2017` masks values absent from the published Li 2017 table and saves
numbered backups. It processes drugs whose fitting outputs exist. hERG binding
and Hill outputs are linked by row number, using sample 0 for optimal inputs
and one-based sample indices for uncertainty inputs. Bootstrap `.npy` arrays
internally use zero-based cell indices and have an accompanying experiments CSV.

For new compounds, use `bootstrap/herg-fit --data PATH`,
`hill-fit -f PATH`, and `ap --cmaxfile PATH --hergpath PATH --hillpath PATH`.
If the compound has no Cmax entry, doses are interpreted as nM. The standard
table is `cipa/data/ap/CiPA_training_drugs.csv`.

To omit a drug effect, use `ap --dropcurrent ICaL` (or hERG, INaL, INa, IKs,
Ito, IK1). `ap --omit` runs all seven omissions. Combine and analyze each
`results/ap/fixed/drop_CURRENT` or `uncertainty/drop_CURRENT` directory
separately. `--no-opt-scaling` disables the Dutta optimized scaling factors.

After simulating the training drugs:

```powershell
.venv\Scripts\python.exe -m cipa tdp-error results/ap/fixed
.venv\Scripts\python.exe -m cipa tdp-error results/ap/uncertainty
```

Classification uses the CiPA labels in the therapeutic-concentration table,
requires all three categories in each training split, and excludes simulations
whose selected beat has missing `max_dv`, following upstream. Output includes
training and cross-validation mean absolute category errors, coefficients,
and per-observation probabilities. Fits that fail to converge issue warnings.

## Python API

```python
from cipa.ap import simulate

metrics, trace = simulate(beats=2, nsave=1)
print(metrics["APD90"], metrics["qNet"])
```

Time is in milliseconds; concentrations are nM; voltage is mV; intracellular
calcium is mM; currents are microA/microF. The `qNet` ODE state and metric CSV
retain upstream's current-times-ms units. Divide qNet by 1,000 for
microC/microF, as done in the qNet plot.

## Conversion map

| Original | Python replacement |
|---|---|
| `hERG_fitting/models/hergmod.c` | `cipa/_herg_equations.py` |
| `AP_simulation/models/newordherg_qNet.c` | `cipa/_ord_equations.py` |
| `generate_bootstrap_samples.R` | `cipa bootstrap` |
| `hERG_fitting.R`, setup and helper functions | `cipa herg-fit`, `cipa/herg.py` |
| `process_boot_results.R` | `cipa herg-summary` |
| `Milnes_sensitivity.R` | `cipa herg-sensitivity` |
| `IC50_mcmc.R` | `cipa hill-fit`, `cipa/hill.py` |
| `compute_mcmc_CI.R` | `cipa hill-summary` |
| `keep_only_Li2017_IC50.R` | `cipa filter-li2017` |
| `process_patchclamp_data.R` | `cipa process-patchclamp` |
| `AP_simulation.R`, parameter and metric helpers | `cipa ap`, `cipa/ap.py`, `cipa/models.py` |
| `combine_results.R` | `cipa combine` |
| `compute_qNet_CI.R` | `cipa qnet-ci` |
| `compute_TdP_error.R` | `cipa tdp-error` |
| Shell launchers | CLI sample/dose ranges, usable on Windows/Linux/macOS |

## Numerical methods and intentional differences

* Equation kernels are generated reproducibly by `tools/translate_models.py`.
  Python expressions retain the original C parameter and state order.
  `exp(n*log(D))` is expressed as `D**n` to handle zero drug concentration.
  The underlying ORd equations retain their upstream removable singularity
  at voltage exactly zero; ordinary solver trajectories cross it without
  evaluating exactly there. Direct evaluations at exactly zero need care.
* AP integration uses SciPy LSODA with upstream tolerances (1e-6), splitting
  each beat at the 0.5 ms stimulus boundary. qNet is reset after every beat.
  The beat with maximum eligible repolarization dv is selected from the saved
  beats; the saved trace and starting state belong to that selected beat.
* At each fixed hERG voltage and concentration the nine-state system is linear.
  The port uses matrix exponentials instead of upstream's loose-tolerance LSODA
  integration. Voltage timings, 180 s preconditioning, and ten 25 s sweeps
  are preserved. The objective retains fractional-current and trapping terms;
  negative-state penalties are accumulated over all concentrations (upstream
  inadvertently kept only the last concentration's penalty).
* hERG fitting uses Python `cma` in the original logarithmic 0-to-10 encoding,
  including the signed Vhalf ratio, with population 80 and tolerance 0.001.
  The optimizer and NumPy RNG are not bit-for-bit equivalents of R's versions.
* Hill fitting uses bounded SciPy least squares with multiple starts.
  MCMC uses adaptive Gaussian Metropolis, uniform bounded parameter priors,
  and inverse-gamma conditional variance updates with prior strength N and
  the fitted residual variance. Adaptation ends at burn-in. It does not
  implement FME's delayed rejection or automatic unbounded burn-in retries.
  Diagnostics use a batch-means early/late Geweke-style z-score; these are
  screening diagnostics, not proof of convergence. Posterior samples and
  fitted values can differ from R's FME implementation.
* Confidence intervals honor the requested alpha (some upstream scripts
  hardcoded 0.025/0.975). Missing bootstrap rows retain their sample positions.
  Aggregation and qNet interval analysis check missing/duplicate samples.
* Patch-clamp preprocessing fixes upstream's `Units`/`units` capitalization
  bug so micromolar values actually convert to nM. The bundled processed
  `drug_block.csv` is preserved exactly; the processing command writes a
  separate output by default.
* TdP regression uses statsmodels proportional-odds logistic regression.
  Its parameterization differs from R rms; exported intercepts and slope
  are converted to the upstream `P(Y >= k) = logistic(intercept_k+slope*qNet)`
  convention. Prediction outputs contain individual probabilities rather
  than upstream's aggregated counts.
* RDS outputs are replaced by CSV, NPY, and JSON. Figures use Matplotlib;
  presentation layouts are not replicas of the R figures. Custom cycle
  lengths are supported but use the supplied CL2000 initial state; pace
  sufficiently before interpreting those runs.

## Validation and provenance

`tests/` covers probability conservation, matrix-exponential versus LSODA
solutions, a physiological control AP, qNet/current integration, known Hill
parameters, reproducible sampling, bootstrap strata, drug scaling/omission,
unit conversion, and missing-sample handling.

`tools/validate_c_models.py` compiles the pinned original C equations with
R interface/timeout code removed and compares their derivatives and current
outputs with the Python port over 60 cases per model, including zero drug,
both stimulus phases, and all three AP cell types. It requires the optional
development compiler `ziglang`. See `validation_c_models.json` for measured
errors. This is equation parity, not proof of identical R statistical results.

`tools/smoke_workflows.py` runs every major workflow with tiny fits and short
AP simulations, writing `results/smoke/`. Those parameter estimates are
deliberately unconverged and must not be used as scientific results.
The classification smoke uses synthetic labelled observations. See
`validation_workflows.json` for the run status. Full 2,000-sample studies
have not been rerun as part of this conversion.

The FDA source authors are Kelly Chang and Zhihua Li; the ORd model credits
Thomas O'Hara, Yoram Rudy, and Washington University in St. Louis. This is a
modified derivative, not an FDA release. Original licenses and disclaimers
are retained in `cipa/licenses/` and the reference checkout. Dataset source
descriptions remain in the reference repository's data README files.
