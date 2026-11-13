BETA-ADRENERGIC SINGLE-CELL STARTER -- EXPLORATORY

Derived from the user-supplied Gong et al. MATLAB source archive.
Paper: Gong et al., J Mol Cell Cardiol 143 (2020) 96-106.
DOI: 10.1016/j.yjmcc.2020.04.009
Source: https://github.com/JQXGong/Ohara-beta-adrenergic
Original authors: Jingqi Q.X. Gong, Monica E. Susilo, Anna Sher,
Cynthia J. Musante, Eric A. Sobie. Underlying models: O'Hara et al.
and Heijman et al. Original source files retained in author_source/.
Python translation and starter runner prepared with ChatGPT/Codex.
Acknowledge and check this assistance as required by your course.

Requirements: Python, NumPy, SciPy; Matplotlib for plotting.
Import: from gong_beta import run_cell
result = run_cell(iso=0.0, beats=5, pcl=1000)

This translates the AUTHOR'S integrated model. It does not modify,
replace, or integrate your teammate's CiPA Markov-IKr model.
Author conductance scaling is enabled; epicardial cell type is used.
58 electrical states + 57 signaling states = 115 states.
Input ISO is micromolar; time is ms; voltage mV; internal calcium mM.
Output calcium_uM is converted to micromolar for plotting.

Signaling is independently preconditioned for 10,000 seconds at the chosen
ISO before electrical pacing. Full signaling/electrical equations are then
integrated together with BDF, splitting at stimulus boundaries.
Stimulus: -80 uA/uF for 0.5 ms, starting 10 ms into each beat.
First electrical states are the author's unconditioned defaults.
Five-beat runs are operational checks, not steady-state validation.
Return 'state' can initialize a later continuation run via state=... .

APD90 in this helper is interpolated at 90% repolarization relative to
pre-stimulus voltage and measured from stimulus onset. It is not guaranteed
to match the author's findAPD convention. Multiple crossings/EADs require
separate analysis; this helper reports the first downward crossing.

Checks performed: finite 115-state derivatives at ISO 0 and 1; successful
five-beat integrated runs; finite states and positive Na/K/Ca concentrations;
one-beat tolerance sensitivity check. validation_results.json records short-run
results. These checks DO NOT establish MATLAB equation/output parity,
steady state, drug compatibility, or reproduction of published figures.

Translation follows the supplied equations and uses an isolated namespace
per cell model. Source conditional square roots use NumPy scimath in the
signaling routines to mirror MATLAB evaluating unused complex branches.
No deliberate physiological simplifications or calibration were introduced.
Do not use this preliminary implementation for final project conclusions
until it is compared against original MATLAB and conditioned sufficiently.
