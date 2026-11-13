"""124-state hybrid using the existing CiPA cell interface."""
import numpy as np
from scipy.integrate import solve_ivp
from . import signaling, _equations
from cipa.cell import ORdCell
from cipa.models import AP_STATES, CURRENTS, initial_ap, parameter_vector

EXTRA_STATES = ('hPf', 'hPs', 'jP', 'hBPf', 'hBPs', 'jBP', 'dP', 'fPf', 'fPs',
                'fcaPf', 'fcaPs', 'fBPf', 'fcaBPf', 'xs1P', 'xs2P', 'JrelP', 'JrelBP')
CONTROL_NAMES = ('fICaLP', 'fIKsP', 'fPLBP', 'fTnIP', 'fINaP', 'fINaKP', 'fRyRP', 'fIKurP')
SIGNAL_OFFSET = len(AP_STATES) + len(EXTRA_STATES)


class CiPABetaCell(ORdCell):
    name = 'cipa-gong-beta'
    state_names = AP_STATES + EXTRA_STATES + tuple('beta_'+n for n in signaling.STATE_NAMES)
    output_names = CURRENTS + CONTROL_NAMES + ('Whole_cell_PP1',)

    def __init__(self, iso=0., regulation=True, precondition_seconds=10000.):
        if not np.isfinite(iso) or iso < 0:
            raise ValueError('ISO must be finite, nonnegative, and specified in micromolar')
        if not np.isfinite(precondition_seconds) or precondition_seconds < 0:
            raise ValueError('precondition_seconds must be finite and nonnegative')
        self.iso = float(iso)
        self.regulation = bool(regulation)
        self.precondition_seconds = float(precondition_seconds)
        self.constants = signaling.constants(self.iso)
        self._signaling_initial = None

    def signal_initial_state(self):
        if self._signaling_initial is None:
            state = signaling.initial_state()
            if self.regulation and self.precondition_seconds:
                sol = solve_ivp(lambda t, y: signaling.signaling(t, y, self.constants),
                    (0., self.precondition_seconds*1000.), state, method='BDF', rtol=1e-6, atol=1e-9)
                if not sol.success or not np.isfinite(sol.y).all():
                    raise RuntimeError(f'Signaling preconditioning failed: {sol.message}')
                state = sol.y[:, -1]
            self._signaling_initial = state.copy()
        return self._signaling_initial.copy()

    def initial_state(self, concentration=0):
        base = initial_ap(concentration)
        # Initialize additional gates from their instantaneous steady states,
        # and release populations consistently with CiPA's starting release.
        v = base[0]
        hp = 1./(1.+np.exp((v+82.90+5.)/6.086))
        hbp = 1./(1.+np.exp((v+89.1+5.)/6.086))
        dp = 1./(1.+np.exp(-(v+3.940+12.)/4.230))
        fp = 1./(1.+np.exp((v+19.58+8.)/3.696))
        xs = 1./(1.+np.exp(-(v+11.60)/8.932))
        extra = np.array([hp, hp, hp, hbp, hbp, hbp, dp, fp, fp, fp, fp, fp, fp,
                          xs, xs, 2.5*base[36], 2.5*base[37]])
        return np.r_[base, extra, self.signal_initial_state()]

    def controls(self, state):
        if not self.regulation:
            return np.zeros(8), .1371
        settings = signaling.effective(state, signaling.SimpleNamespace(const_signaling=self.constants))
        c = self.constants
        pp_sum = c.pp1_K-c.PP1_cyt+state[38]
        free = .5*(np.sqrt(pp_sum**2+4*c.pp1_K*c.PP1_cyt)-pp_sum)
        pp1 = c.PP1_cav/c.vr_cav+c.PP1_eca/c.vr_eca+free/c.vr_cyt
        return np.array([getattr(settings, n) for n in CONTROL_NAMES]), pp1

    def evaluate(self, time, state, parameters, stimulating):
        p = parameter_vector(parameters) if isinstance(parameters, dict) else parameters
        beta = state[SIGNAL_OFFSET:]
        control, pp1 = self.controls(beta)
        electrical, currents = _equations.evaluate(0. if stimulating else 1., state[:SIGNAL_OFFSET],
                                                   p, control, pp1, self.regulation)
        signal = signaling.signaling(time, beta, self.constants) if self.regulation else np.zeros(57)
        return np.r_[electrical, signal], np.r_[currents, control, pp1]

    def metadata(self):
        return dict(beta_adrenergic_model='Gong 2020 signaling and substrate regulation',
            iso_uM=self.iso, regulation_enabled=self.regulation,
            precondition_seconds=self.precondition_seconds, iso_protocol='constant from preconditioning onward',
            initial_signal_max_abs_derivative=float(np.max(np.abs(signaling.signaling(
                0., self.signal_initial_state(), self.constants)))),
            celltype='CiPA parameter celltype (default endocardial)',
            scaling='CiPA; no Gong baseline scaling applied',
            drug_assumption='Same Hill block for all phosphorylation populations; dynamic CiPA hERG retained')
