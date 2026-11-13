# Preparing the quinidine study for beta-adrenergic stimulation

## Scientific scope

Heijman, Volders, Westra and Rudy (2011), *Local control of beta-adrenergic
stimulation: Effects on ventricular myocyte electrophysiology and Ca2+-transient*:

- Paper: https://pubmed.ncbi.nlm.nih.gov/21345340/
- Full text and supplement: https://pmc.ncbi.nlm.nih.gov/articles/PMC3075371/
- Authors' model download page: https://rudylab.wustl.edu/code-downloads/
- Additional model representation: https://models.cellml.org/workspace/6ec

The paper describes canine epicardial electrophysiology coupled to local
caveolar, extracaveolar, and cytosolic signaling domains. It includes beta1/beta2
receptors, cAMP/PKA, phosphorylation-dependent substrate populations, and
CaMKII interactions. The present CiPA equations describe human ventricular
electrophysiology. Adding the entire canine electrophysiology model would
replace the human model rather than simply stimulate it.

The current refactor prepares a human ORd plus signaling hybrid. It does not
implement or claim validation of beta-adrenergic effects. No ISO flag or invented
conductance multipliers are exposed. The downloaded paper's code has not yet
been imported; source licenses, exact version, and baseline reproduction must
be checked as part of that integration.

## Architecture now available

| Layer | File | Responsibility |
|---|---|---|
| Study scope | `cipa/study.py` | Quinidine and control; dedicated Cmax table |
| Drug effects | `cipa/models.py` | Load hERG/Hill fits and calculate channel block |
| Cell interface | `cipa/cell.py` | Named states/outputs, evaluate, initial state, beat reset |
| Baseline equations | `cipa/_ord_equations.py` | Existing generated human CiPA kernel |
| Pacing and integration | `cipa/ap.py` | Stimulus boundaries, absolute time, state persistence |
| Metrics/figures | `cipa/analysis.py` | AP/calcium/qNet analysis and plots |

Pass a CellModel-compatible object using `simulate(..., model=my_model)`.
`simulate_beat(..., model=my_model, start_time=...)` is also available.
The runner supports extra signaling states and additional output columns,
and calls `reset_beat` instead of resetting the last state unconditionally.
Only qNet is reset by ORdCell; future phosphorylation and messenger states
must persist. The runner supplies absolute milliseconds to `evaluate`, even
when replaying the selected beat. The explicit `stimulating` argument carries
the fast pacing pulse independently of the signaling clock.

The default adapter consumes the existing ORd parameter vector. A hybrid can
reuse this vector for baseline parameters and store signaling configuration
in its own constructor. Initial parameters remain drug-adjusted by
`ap_parameters`. The metrics contract currently requires the named columns
`v`, `cai`, `qNet`, and `dv`; new models must supply these. Core state names
must remain unique and output names must not duplicate state names.

## Next implementation steps

1. Obtain and pin the authors' source; preserve license and attribution.
   Reproduce its no-ISO and ISO canine baseline before extracting signaling.
2. Port signaling states, reaction equations, compartment volumes, and
   initial conditions independently, recording concentration and time units.
3. Map each substrate to the human equations: ICaL, IKs, calcium reuptake/
   PLB/SERCA, RyR release, Na/K pumping, and any other supported targets.
   Identify canine substrates absent from ORd rather than adding them silently.
4. Reconcile CaMKII: ORd already has CaMK-dependent regulation. Avoid counting
   the same phosphorylation pathway twice. Define which calcium compartment
   drives each reaction and which baseline phosphorylation fractions are used.
5. Implement an explicit hybrid equation module. Reweight phosphorylated and
   non-phosphorylated populations where required; simply multiplying reported
   currents afterward would leave voltage and concentration derivatives wrong.
   Apply quinidine block consistently to the appropriate modified currents.
6. Add an ISO concentration/onset protocol with explicit units and timing.
   Treat onset/offset as integration boundaries, just like the pacing stimulus.
7. Validate no-ISO hybrid behavior against the current human baseline, then
   check ISO dose responses, phosphorylation, APD, calcium amplitude/relaxation,
   and pacing-rate dependence. Compare control, quinidine, ISO, and both.

The original multi-drug TdP ordinal classifier cannot be trained from quinidine
alone. The active CLI rejects it; compare continuous APD, calcium, qNet, and
afterdepolarization measures within this study instead. Published CiPA risk
thresholds cannot be assumed valid for a new human/signaling hybrid.

## Outputs and reproducibility

Each dose now includes `run_config.json` with model identity, concentration,
pacing, sample, and input paths. New hybrid models should record their source
revision, ISO protocol, and coupling parameters in an accompanying metadata file.
Use separate `--output` roots for baseline and future ISO experiments to avoid
overwriting or aggregating unlike experimental conditions.
