# CiPA + Gong beta-adrenergic single-cell model

This package combines the existing CiPA human ventricular model and quinidine
pharmacology with the partner-supplied `beta-signaling/gong_beta.py` translation
of Gong et al. (2020). The existing `cipa/` and `beta-signaling/` implementations
are unchanged. Fitting still uses `python -m cipa`; combined AP simulations use
`python -m cipa_beta`.

It uses the project's existing virtual environment and dependencies. Run from
the repository root; if you want module discovery outside that directory,
refresh the editable installation with `python -m pip install -e . --no-deps`.

## Run from the project root

First obtain quinidine fits using your existing workflow. Existing fits are
reused, not recomputed or recalibrated by this package:

```powershell
python -m cipa herg-fit -d quinidine
python -m cipa hill-fit -d quinidine
python -m cipa filter-li2017
```

Short operational checks (five beats are not steady-state results):

```powershell
# No drug, basal beta signaling
python -m cipa_beta ap -d control --iso 0 --beats 5 --nsave 2 --trace-beats 5

# No drug, 1 micromolar ISO
python -m cipa_beta ap -d control --iso 1 --beats 5 --nsave 2 --trace-beats 5

# Quinidine at 1x Cmax, with ISO
python -m cipa_beta ap -d quinidine -x 1 --iso 1 --beats 5 --nsave 2 --trace-beats 5
```

Full-length pacing with ten exported cycles:

```powershell
python -m cipa_beta ap -d quinidine -x "1-4" --iso 1 --trace-beats 10
```

Defaults are 1,000 beats, 2,000 ms cycle length, 250 beats analyzed, 1 ms
output sampling, endocardial cell parameters, and optimized CiPA scaling.
ISO is **micromolar**; quinidine Cmax is **3237 nanomolar**. These are different
units. Quinidine dose multipliers, sample indices, `--hergpath`, `--hillpath`,
`--dropcurrent`, and `--no-opt-scaling` follow the existing CiPA workflow.

For uncertainty runs, supply previously fitted hERG bootstrap and Hill MCMC
rows with `-i "1-2000"`. The number of available rows must match. Convergence
warnings from fitting still need to be assessed before interpreting uncertainty.
To use the longer Hill run, specify `--hillpath results/hill_long`; its filtering
can be performed with `python -m cipa filter-li2017 --output results/hill_long`.

## Outputs

By default the two ISO conditions are stored separately:

```text
results/cipa_beta/
  iso_0_uM/fixed/control/...
  iso_1_uM/fixed/control/...
  iso_1_uM/fixed/quinidine/
    metrics.csv
    dose_1/
      trace.csv
      trace.png
      trace_cycles.csv          # when --trace-beats > 0
      trace_cycles.png
      beat_metrics.csv
      states_max_dv.csv
      checkpoints.csv
      run_config.json
```

Uncertainty runs use `uncertainty/DRUG/SAMPLE/`. `--no-regulation` uses the
`regulation_off/` condition. `--output` changes the root. Reusing the same
condition/dose paths overwrites that run, so use distinct roots for different
pacing or initialization experiments.

The existing plot scripts can read these trace CSVs by changing their file
paths. Voltage, calcium and the seven currents retain their existing names.
Extra columns record all 57 signaling states with `beta_` prefixes, the 17
new electrical gates/release states, eight effective phosphorylation controls,
and `Whole_cell_PP1`. `fIKurP` controls IKb in the supplied Gong equations;
the hybrid does not add a new IKur current.

The single trace is selected with the existing CiPA maximum-repolarization-dv
rule; the cycle file contains the final N consecutive beats. qNet is computed
from the same six modified currents and retains the original units; divide
by 1000 for microC/microF. qNet is reset each beat; signaling states are not.

## What was combined

- CiPA's 50 states, geometry, ion concentrations, physical constants, baseline
  conductance scaling, and dynamic Markov hERG drug-binding/trapping equations.
- Gong's 57 signaling equations and normalized effective phosphorylation
  mapping, with local named state indices instead of its 58-state offset.
- Seventeen additional electrical states for phosphorylated INa, ICaL, IKs,
  and RyR populations; 124 total states.
- Gong regulation of INa, ICaL/ICaNa/ICaK, IKs, PLB/SERCA, RyR, TnI buffering,
  Na/K pumping, background IKb, and PP1-dependent CaMK deactivation.
- PKA-only, CaMK-only, doubly phosphorylated, and unphosphorylated populations
  replace the corresponding CiPA expressions. They are not counted twice.
- CiPA Hill conductance/permeability factors act on all populations of each
  affected channel. Equal quinidine sensitivity across phosphorylation states
  is an explicit assumption; this model does not establish it experimentally.

Gong's alternative baseline scaling, epicardial default, Hodgkin-Huxley IKr,
and 10 ms stimulus offset were not imported. The hybrid retains CiPA pacing
at beat time zero. The supplied Gong 0.025 mM Cass ceiling in the ICaL GHK
expression is also excluded to retain the CiPA flux equations and disabled
baseline parity. All other transferred regulation expressions follow the
local partner-supplied implementation.

`tools/build_cipa_beta.py` extracts the signaling and transplants the reviewed
regulation blocks into a separate CiPA kernel. Generated `signaling.py` and
`_equations.py` are ordinary Python modules, with optional Numba for the
electrical kernel. `provenance.json` records source SHA256 hashes. Regenerate
only when intentionally updating the reviewed sources; it overwrites those
generated files.

## Initialization and numerical methods

ISO is constant throughout each run. Following the supplied starter, the 57
signaling states are first integrated alone for **10,000 seconds** at that ISO.
Use `--precondition-seconds` to change this (zero skips preconditioning).
This represents ISO established before electrical pacing, not acute ISO
addition during the recording. Dynamic onset/offset protocols are not included.

CiPA's conditioned CL2000 electrical state initializes the common states.
New gates start at instantaneous steady states at the initial voltage; new
release populations start at 2.5x their CiPA counterparts. Electrical states
must still be paced to settle into the hybrid condition. Increasing the
number of beats does not itself prove steady state.

All 124 states are solved together with BDF, rtol 1e-6 and atol 1e-9, with an
explicit split at the 0.5 ms pacing stimulus boundary. Signaling rates already
include the conversion to milliseconds. They are not converted a second time.
Absolute experiment time persists across beats and selected-beat replay.
Metadata includes ISO, preconditioning duration and its terminal maximum
absolute signaling derivative, pacing, model identity, drug concentration,
and parameter-input paths.

For an explicit regression condition without beta regulation:

```powershell
python -m cipa_beta ap -d control --no-regulation --beats 5 --nsave 2
```

Disabled regulation sets all effective beta controls to zero, freezes the
signaling states and restores CiPA's CaMK deactivation. It retains the CiPA
common-state equations to floating-point precision. **ISO=0 with regulation
enabled is a separate basal-signaling condition**, not an instruction to
disable the signaling model.

## Validation and limits

Run `python -m pytest tests/test_cipa_beta.py -q`.

Checks cover disabled-regulation equation parity across all three cell types
with nonzero drug-binding parameters, hERG probability conservation, active
regulation against the supplied Gong electrophysiology with matched CiPA
constants and IKr, exact signaling extraction, equal Hill block across
phosphorylated populations, finite coupled trajectories, cycle export, and
signaling persistence across beat resets.

The full regression suite passed 21 tests. `tools/validate_cipa_beta.py` also
ran five beats for control and quinidine at ISO=0 and ISO=1 uM, checking finite
states, positive ionic concentrations, hERG probability conservation and cycle
exports. See `validation.json`; its plot is saved under
`results/cipa_beta_validation/comparison.png`. The quinidine+ISO short run
exhibited afterdepolarizations, and ISO runs produced negative qNet. These
observations are exploratory model outputs, not validated risk predictions.

These are implementation checks. The partner's translation does not yet claim
MATLAB parity or reproduction of published figures. The human hybrid has not
been recalibrated or validated as a CiPA risk classifier; published qNet risk
thresholds must not be assumed transferable. Basal/ISO/drug/both conditions
should be compared after adequate conditioning and further numerical and
physiological validation. No phosphorylation-dependent drug affinity or
ISO-induced heart-rate change is assumed.

See `NOTICE`, `provenance.json`, `../cipa/licenses/`, and the partner's
`../beta-signaling/README.txt` for attribution and source limitations.
