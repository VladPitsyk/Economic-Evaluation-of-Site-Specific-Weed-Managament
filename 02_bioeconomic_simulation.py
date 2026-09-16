"""
Bio-economic simulation of site-specific weed management under imperfect sensing.


Description:
    Evaluates the annualised net return of four site-specific herbicide
    application strategies relative to a uniform broadcast application, for a
    single field at a single spraying resolution. The model operates on a
    rotated grid layer carrying per-cell weed detection counts (produced by
    01_grid_creation_weed_counts.py) and runs inside the QGIS Python console,
    which supplies the qgis.core API used to read and write the layer.

    For every grid cell the model:
      1. corrects the detected weed counts for the sensor recall and allocates
         the undetected weeds between detected-weed patches and otherwise
         empty cells (Eq. 7);
      2. applies the decision rule of the selected strategy to obtain a spray
         decision and dose;
      3. computes the Year-0 yield loss from surviving weeds using the
         rectangular-hyperbola competition model (Eqs. 1-5);
      4. projects the seedbank over a 10-year horizon under probabilistic,
         density-dependent detection in each future year, and discounts the
         resulting herbicide and yield-loss stream to a future weed penalty
         (Eqs. 8-12);
      5. computes the marginal shadow cost of one additional weed of each
         class, holding the cell's spray decision fixed.

    Results are appended to a CSV, one row per field x resolution x strategy,
    together with the sensitivity analyses reported in the paper and the run
    configuration that produced them. Per-cell values are written back to the
    grid layer for the marginal shadow cost figure.

Inputs:
    - A QGIS project containing a grid layer (default name 'grid_rotated')
      with the fields DicotCount and MonoCount, as written by
      01_grid_creation_weed_counts.py.

Outputs:
    - <RESULTS_CSV_PATH>                    main results, one row per run
    - <RESULTS_CSV_PATH>_recall_sweep.csv   net return by sensor recall
    - <RESULTS_CSV_PATH>_dd_kill_sweep.csv  net return by background-dose efficacy
    - Per-cell attributes written to the grid layer, including Marg_SC_D and
      Marg_SC_M, which are exported to CSV as input to
      04_marginal_shadow_costs.R.

Requirements:
    QGIS >= 3.28 with Python >= 3.9 (run from the QGIS Python console).
    No third-party packages beyond those bundled with QGIS.

Usage:
    1. Open the QGIS project for one field at one resolution.
    2. Set RESULTS_CSV_PATH and SCENARIO_ID in the USER CONFIGURATION section.
    3. Run the script from the QGIS Python console.
    4. Repeat for each field, resolution and strategy.

    All CSV outputs are opened in APPEND mode, so a full set of results is
    accumulated across runs. Archive or delete the CSVs before re-running a
    strategy, otherwise old and new rows will coexist; the configuration
    columns written with every row identify which settings produced it.

Notes on reproducing the published results:
    The dual-dose strategy is reported at a background-dose kill fraction of
    DD_BASE_KILL = 0.85 with DD_FUTURE_COST_MODE = 'dose_scaled'. The
    alternative setting 'full_rate' charges the full herbicide cost in
    background-dose years; combined with DD_BASE_KILL = 0.95 it makes the
    dual-dose strategy algebraically identical to broadcast application, which
    is retained only as a validation identity and must not be used to produce
    reported results.
"""

import os
import csv
import math
from qgis.core import QgsProject, QgsField
from qgis.PyQt.QtCore import QVariant
from qgis.utils import iface

# ================= USER CONFIGURATION =================

SCENARIO_ID = 1
scenario_lookup = {
    1: 'SCENARIO_1_ECON',    # ET  — myopic economic threshold (Year 0 only)
    2: 'SCENARIO_2_SPOT',    # PB  — presence based
    3: 'SCENARIO_3_SECURE',  # DD  — dual dose
    4: 'SCENARIO_4_DYNET',   # forward-looking ET
}
SCENARIO = scenario_lookup[SCENARIO_ID]

TECH_SCENARIO = 'BASED'

# Absolute path to the main results CSV. The two sweep files are written
# alongside it with '_recall_sweep' and '_dd_kill_sweep' appended.
RESULTS_CSV_PATH = r'C:\path\to\output\results_thesis.csv'

# ----- TOGGLES -----
USE_CELL_CACHE = True    # leave True
ENABLE_RECALL_SWEEP = True    # net return recomputed at each recall value
ENABLE_DD_KILL_SWEEP = True    # runs automatically on DD runs only
ANNUALISE_SHADOW_COST = False   # report marginal shadow costs as 10-year NPV

# DUAL-DOSE FUTURE CHEMICAL COST — affects the DD strategy only.
#   'dose_scaled' (DEFAULT, used for all reported results) charges
#                 MAX_HERB_COST_YR * DD_BASE_DOSE
#                 in background-dose years and the full cost in full-dose
#                 years. This is what the DD strategy as described in the
#                 manuscript actually does, and is consistent with how ET/PB
#                 future costs are handled (cost follows the applied dose).
#   'full_rate'   charges DD the FULL herbicide cost in every future year, even
#                 in background-dose cells. Combined with DD_BASE_KILL = 0.95
#                 this makes DD algebraically identical to broadcast (the
#                 validation identity: DD differential = 0 everywhere).
#                 VALIDATION USE ONLY — do not report results from 'full_rate'.
DD_FUTURE_COST_MODE = 'dose_scaled'

# ET DECISION RULE IN FUTURE YEARS
#
# In Year 0 the economic threshold is evaluated on the INFERRED TRUE density
# (bio_dens = T/area), i.e. detections scaled up by the known false-negative
# rate. An earlier formulation evaluated the same threshold in years 1-10 on the
# RAW DETECTED density (emerged * recall) with no scale-up, so the strategy used
# different information in year 0 than in later years.
#
#   'full_rate'        keeps the original behaviour (raw detected). Use for validation.
#   'consistent' applies the same scale-up factor used in Year 0, so the
#                economic threshold is evaluated on the estimated true density
#                in every year.
#
# Affects ET (SCENARIO_1_ECON) and DynET (SCENARIO_4_DYNET) ONLY.
# PB and DD are untouched and must stay bit-identical — that is the validation.
ET_FUTURE_RULE = 'consistent'

# =============================================================================
# WHAT INFORMATION THE ECONOMIC THRESHOLD ACTS ON
# =============================================================================
#
#   'detected'  The threshold is applied to the weed counts the sensor actually
#               reports. This represents a commercial system as an operator
#               would use it: the recall of the detection algorithm is not known
#               to the farmer, is not measured in the field, and differs between
#               systems, so it cannot enter the decision rule.
#
#   'true'      The threshold is applied to detected counts scaled up by the
#               known false-negative rate, i.e. to the estimated TRUE density.
#               This represents an operator who knows the recall of their system
#               and corrects for it.
#
# Applies to SCENARIO_1_ECON only, in Year 0 AND in years 1-10.
#
# SCENARIO_4_DYNET is deliberately NOT affected: it is defined as the informed
# strategy, using the estimated true density and the discounted future weed
# penalty. The two together give a ladder of increasing information:
#
#     PB     treat any cell with a detection            (no threshold)
#     ET     economic threshold on what the sensor sees (current-season only)
#     DynET  economic threshold on estimated true       (current + future)
#            density, including the avoided future
#            weed penalty
#
# NOTE ON RESOLUTION: at 50 x 50 cm and 25 x 25 cm a single detected weed
# already clears the threshold under either setting, so ET still converges on PB
# at those resolutions. Only the 1 x 1 m case changes, where the threshold moves
# from 2 detected weeds to 3.
ET_INFORMATION = 'detected'

# =============================================================================
# DETECTION IN FUTURE YEARS  (years 1-10)
# =============================================================================
#
#   'expected'      Detected weeds are computed as
#                   emerged x recall, a continuous expected value that is above
#                   zero whenever ANY weed emerges. Every cell containing weeds
#                   is therefore identified and treated from year 1 onward, so
#                   site-specific strategies behave like broadcast application
#                   after Year 0 and escapes never compound. Use for validation.
#
#   'probabilistic' Detection is modelled as a random process, so cells holding
#                   few weeds are repeatedly missed and their seedbank
#                   compounds, as it would in the field. See the block comment
#                   above _future_liability() for the full derivation.
#
# NOTE: Year 0 is NOT affected by this setting. The Year-0 spray map comes from
# the actual UAV detections (P_i), which are observed data and already embody
# the real detector's performance. Only years 1-10, which are simulated, need a
# detection model.
FUTURE_DETECTION = 'probabilistic'

# Branching control (used only when FUTURE_DETECTION = 'probabilistic').
# Each year a cell either is or is not treated, so a 10-year horizon has up to
# 2^10 = 1024 trajectories. Paths carrying less probability than PRUNE_PROB are
# dropped, and trajectories that have converged to near-identical seedbanks are
# merged, which keeps the state count small without materially affecting the
# expected value. Diagnostics are printed at the end of the run.
# Convergence was verified: the expected liability was identical to four
# decimal places for MAX_STATES between 128 and 4096 and MERGE_DIGITS between 3
# and 5, so the pruned trajectories carry negligible probability. The defaults
# below sit comfortably above the natural state count (~600).
PRUNE_PROB = 1e-6
MERGE_DIGITS = 3      # significant digits used when merging similar seedbanks
MAX_STATES = 1024   # ceiling on simultaneous trajectories

# Fraction of undetected weeds allocated to cells that already contain
# detections (the remainder is spread over zero-detection cells). Used both for
# the Year-0 true-count estimate and, under ET_FUTURE_RULE = 'consistent', for
# the future-year decision scale-up.
PATCH_FRACTION = 0.70

tech_cost_lookup = {
    'BASED':      0,
    'SCENARIO 1': 0,
    'SCENARIO 2': 0,
    'SCENARIO 3': 0,
    'SCENARIO 4': 0,
}
sswm_annuity_per_ha = tech_cost_lookup[TECH_SCENARIO]

GRID_LAYER_NAME = 'grid_rotated'
INPUT_DICOT_FIELD = 'DicotCount'
INPUT_MONOCOT_FIELD = 'MonoCount'

# ----- CURRENT YEAR PARAMETERS -----
dicot_competition_index = 8.0
monocot_competition_index = 2.0
weed_height_modifier = 0.75

yield_loss_initial_slope = 0.50
yield_loss_asymptote_pct = 60.0

weed_free_yield_kg_per_m2 = 1.015
crop_price_eur_per_kg = 0.229
herbicide_cost_fullrate_eur_per_m2 = 0.0130
full_kill_rate = 0.95

DD_BASE_DOSE = 0.5
# Kill fraction of the reduced BACKGROUND dose (full-dose cells always
# use full_kill_rate = 0.95). 0.85 is the manuscript's main-results assumption
# (mid-range of published reduced-rate efficacies); it applies both to the
# Year-0 background cells and, via calculate_10yr_liability_v2, to
# background-dose branches in years 1-10. Set 0.95 only to reproduce the
# DD == broadcast validation identity together with DD_FUTURE_COST_MODE='full_rate'.
DD_BASE_KILL = 0.85
security_base_dose = DD_BASE_DOSE
security_base_kill = DD_BASE_KILL

# ----- DETECTOR QUALITY -----
recall_dicot = 0.41
recall_monocot = 0.41

# ================= BIOLOGICAL ENGINE =================
DISCOUNT_RATE = 0.04
SIM_YEARS = 10
MAX_HERB_COST_YR = 130.0
# SEED SET OF HERBICIDE SURVIVORS (gamma in Eq. 9)
#
# CALIBRATED, not taken directly from the literature. Reported values span
# roughly 0.05-0.40 depending on growth stage at application, mode of action and
# species, so the literature constrains the range but does not fix the value.
#
# 0.13 was chosen because it is the value at which the model reproduces the
# system it represents: run forward under continuous broadcast application at
# 95% efficacy - the management these fields have actually received - the weed
# population settles at 17.3 plants m-2, against an observed median of
# 16.4 plants m-2 across the 19 fields. At the previous value of 0.33 the same
# projection ran away to 119 plants m-2, i.e. the model had no stable state even
# under effective conventional control.
#
# The model is close to bistable in this region: 0.16 gives 33 plants m-2, 0.11
# gives 7, and below about 0.10 the population is driven to extinction. This
# sensitivity should be stated in the manuscript.
#
# NOTE: this parameter governs SEED PRODUCTION only. The separate reduction in
# the COMPETITIVE ability of a stunted survivor is the hardcoded factor of 1/3
# applied to TCL_surv further down, and is unchanged. Equation 5 and Equation 9
# of the manuscript currently use the same symbol (gamma) for both; they need to
# be distinguished now that they take different values.
STUNTED_SEED_PENALTY = 0.13

# derived; equals 8.1109 for r = 0.04, T = 10
ANNUITY_FACTOR = 8.11   # PINNED to match previously computed results.
# Derived value would be (1-(1+r)^-T)/r = 8.110896; pinning avoids 0.011%
# differences between re-run and non-re-run numbers.

SPLIT_SCENARIOS = {
    'split_90_10': 0.90,
    'split_70_30': 0.70,
    'split_50_50': 0.50,
}

# Rebased around the calibrated baseline of 0.13, keeping the same
# relative spread as the original 0.20 / 0.33 / 0.40 (x0.61 and x1.21).
# WITHOUT this rebasing every efficacy scenario would still run away: at the old
# values the equilibrium densities under broadcast are 85, 119 and 0 plants m-2.
KILL_SCENARIOS = {
    'k_low':  {'kill': 0.85, 'stunt': 0.08},
    'k_base': {'kill': 0.95, 'stunt': 0.13},
    'k_high': {'kill': 0.99, 'stunt': 0.16},
}

RECALL_SCENARIOS = {
    'recall_low':  round(recall_dicot * 0.8, 4),
    'recall_base': recall_dicot,
    'recall_high': round(recall_dicot * 1.2, 4),
}

# ----- RECALL SWEEP -----
# The Year-0 spray map is held FIXED at the observed detections, so this
# isolates "fewer weeds escape within cells you already spray" and EXCLUDES
# "you would also find and treat patches you currently miss entirely".
# The break-even recall obtained here is therefore an UPPER bound on the recall
# actually required. Say so explicitly in the paper — it turns a simplification
# into a conservative claim.
RECALL_SWEEP = [0.328, 0.41, 0.492]

# ----- DUAL-DOSE BACKGROUND KILL SWEEP -----
# k_background acts through TWO channels: (a) YEAR 0 — undetected weeds in
# background-dose cells survive at (1 - k_background) instead of 0.05, causing
# extra Year-0 yield loss and extra seed return into the bank; and (b) YEARS
# 1-10 — the sweep value is passed as dd_base_kill into
# calculate_10yr_liability_v2, so background-dose branches under probabilistic
# future detection also kill at k_background. Feeds Fig. 6 / Section 3.5.
# Under 'dose_scaled' the k_bg = 0.95 point is no longer zero by construction:
# it retains the future chemical saving of the background dose.
DD_KILL_SWEEP = [0.95, 0.90, 0.85, 0.80, 0.70]

# ----- WEED BIOLOGY PARAMETERS -----
# 'yl_slope' / 'yl_cap' ARE NEVER USED — kept only so the dicts stay comparable
# Yield loss is computed from the GLOBAL yield_loss_initial_slope
# (0.50) and yield_loss_asymptote_pct (60.0) applied to the COMBINED two-class
# TCL. A single combined competitive load can only carry one asymptote, so this
# is a defensible design — but Table 3 of the manuscript currently describes
# parameters that do not enter the model, including the 35% monocot asymptote.
# ACTION: fix the manuscript text, not the code. The effective per-class slope
# is 0.5 * CI * CIM = 3.00 (dicot) and 0.75 (monocot), and the +/-20%
# sensitivity varies FIVE parameters, not seven.
#
# 'crop_red' is the OVERWINTER SURVIVAL coefficient (phi in Eq. 9), not a crop
# reduction term. Worth renaming before publishing the code as supplementary.
p_dicot = {
    'name':          'dicot_base',
    'max_seeds':     50000,
    'emergence':     0.05,
    'bank_survival': 0.40,
    'intra_comp':    0.41,
    'crop_red':      0.05,
    'yl_slope':      4.0,            # UNUSED
    'yl_cap':        60.0            # UNUSED
}

p_mono = {
    'name':          'mono_base',
    'max_seeds':     2500,
    'emergence':     0.20,
    'bank_survival': 0.50,
    'intra_comp':    0.25,
    'crop_red':      0.10,
    'yl_slope':      1.0,            # UNUSED
    'yl_cap':        35.0            # UNUSED
}


def get_sensitivity_params(p_base, multiplier, tag):
    return {
        'name':          p_base['name'] + '_' + tag,
        'max_seeds':     p_base['max_seeds'] * multiplier,
        'emergence':     min(1.0, p_base['emergence'] * multiplier),
        'bank_survival': min(1.0, p_base['bank_survival'] * multiplier),
        'intra_comp':    p_base['intra_comp'] / multiplier,
        'crop_red':      min(1.0, p_base['crop_red'] * multiplier),
        'yl_slope':      p_base['yl_slope'] * multiplier,
        'yl_cap':        min(100.0, p_base['yl_cap'] * multiplier)
    }


p_dicot_low = get_sensitivity_params(p_dicot, 0.8, 'low')
p_dicot_high = get_sensitivity_params(p_dicot, 1.2, 'high')
p_mono_low = get_sensitivity_params(p_mono,  0.8, 'low')
p_mono_high = get_sensitivity_params(p_mono,  1.2, 'high')


# =====================================================================
# Decision scale-up used by the ET / DynET threshold.
#
# Returns the factor by which a raw detected count is multiplied to obtain the
# estimated true count, using the same patch-allocation logic as Year 0:
#     T / P = 1 + ((1/R) - 1) * PATCH_FRACTION
# At R = 0.41 and PATCH_FRACTION = 0.70 this equals 2.007.
# Returns 1.0 under ET_FUTURE_RULE = 'full_rate', reproducing the original behaviour.
# =====================================================================
def _advance_one_year(bank_d, bank_m, emerged_d, emerged_m, treated,
                      kill_y, seed_penalty_if_treated, chem_cost, p_d, p_m):
    """One year of the population model for a single trajectory.

    Returns (new_bank_d, new_bank_m, cost) where cost is the undiscounted sum of
    herbicide and yield-loss cost for that year. It is
    factored out so the treated and untreated continuations share it exactly.
    """
    if treated:
        survivors_d = emerged_d * (1.0 - kill_y)
        survivors_m = emerged_m * (1.0 - kill_y)
        seed_pen = seed_penalty_if_treated
    else:
        survivors_d = emerged_d
        survivors_m = emerged_m
        seed_pen = 1.0

    TCL_surv = weed_height_modifier * (
        survivors_d * dicot_competition_index +
        survivors_m * monocot_competition_index)
    if treated:
        TCL_surv *= (1.0 / 3.0)
    YL_pct = multispecies_yield_loss_pct(
        TCL_surv, yield_loss_initial_slope, yield_loss_asymptote_pct)
    yl_cost = (YL_pct / 100.0) * (weed_free_yield_kg_per_m2 *
                                  10000) * crop_price_eur_per_kg

    # Seedbank update: bank*survival + new seeds  (ADDITION; Eq. 10 in the
    # manuscript is typeset as a multiplication and needs correcting there)
    seeds_d = (((p_d['max_seeds'] * survivors_d) /
                (1.0 + p_d['intra_comp'] * survivors_d)) * p_d['crop_red']
               if survivors_d > 0 else 0.0) * seed_pen
    seeds_m = (((p_m['max_seeds'] * survivors_m) /
                (1.0 + p_m['intra_comp'] * survivors_m)) * p_m['crop_red']
               if survivors_m > 0 else 0.0) * seed_pen

    return (bank_d * p_d['bank_survival'] + seeds_d,
            bank_m * p_m['bank_survival'] + seeds_m,
            chem_cost + yl_cost)


_BRANCH_STATS = {'max_states': 0}


def _poisson_sf(k, mu):
    """P(X > k) for X ~ Poisson(mu). Used for the detection probability.

    Written out rather than imported so the script has no dependency beyond the
    Python standard library and whatever ships inside QGIS.
    """
    if mu <= 0.0:
        return 0.0
    if k < 0:
        return 1.0
    if mu > 700.0:            # exp(-mu) underflows; detection is certain
        return 1.0
    term = math.exp(-mu)
    cdf = term
    for i in range(1, k + 1):
        term *= mu / i
        cdf += term
        if cdf >= 1.0:
            return 0.0
    return max(0.0, 1.0 - cdf)


def _et_threshold_count(lam_d, lam_m, area, su, scenario='SCENARIO_1_ECON'):
    """Smallest number of DETECTED weeds in a cell that clears the economic
    threshold, given the current class mix.

    Inverts the rectangular hyperbola analytically rather than searching:

        YL* = 100 * herbicide_cost / (yield_goal * price)          target loss
        TCL* = YL* / (slope * (1 - YL*/A))                         from Eq. 3
        x*  = TCL* * area / (CIM * su * mean_competitive_index)     from Eq. 1-2

    Returns None when the threshold cannot be reached at any density, i.e. when
    the herbicide cost exceeds the value of the asymptotic yield loss.
    """
    lam_tot = lam_d + lam_m
    if lam_tot <= 0:
        return None

    value_per_ha = weed_free_yield_kg_per_m2 * 10000.0 * crop_price_eur_per_kg

    # [FIX] DynET counts the discounted carry-over benefit of preventing this
    # cohort from setting seed, not just the current season's loss, so it clears
    # its threshold at a LOWER weed density than the myopic ET. Without this the
    # two strategies are identical in probabilistic mode, because both would
    # fall through to the same threshold.
    effective_cost = MAX_HERB_COST_YR
    if scenario == 'SCENARIO_4_DYNET':
        effective_cost = MAX_HERB_COST_YR / (1.0 + 1.0 / (1.0 + DISCOUNT_RATE))

    yl_target = 100.0 * effective_cost / value_per_ha
    if yl_target >= yield_loss_asymptote_pct:
        return None                      # unreachable: never economic to spray

    denom = yield_loss_initial_slope * \
        (1.0 - yl_target / yield_loss_asymptote_pct)
    if denom <= 0:
        return None
    tcl_target = yl_target / denom

    # competitive index of the average weed in this cell, given the class mix
    frac_d = lam_d / lam_tot
    ci_mean = (frac_d * dicot_competition_index +
               (1.0 - frac_d) * monocot_competition_index)
    if ci_mean <= 0:
        return None

    x_star = tcl_target * area / (weed_height_modifier * su * ci_mean)
    return max(1, int(math.ceil(x_star)))


def _p_treated(lam_d, lam_m, area, recall, scenario, su):
    """Probability that a cell is treated in a given future year.

    THE DETECTION MODEL
    -------------------
    The number of weeds standing in a cell is not known exactly; the simulation
    carries an expected density. Weeds are taken to be distributed at random
    within the cell, so the count N follows a Poisson distribution with mean
    lambda = density x cell area. The detector finds each weed independently
    with probability R (the empirical recall, 0.41). By the thinning property of
    the Poisson distribution the number DETECTED is then exactly

        X ~ Poisson(lambda * R)

    which handles fractional expected counts naturally, unlike a binomial.

    A presence-based sprayer treats the cell if it sees at least one weed:

        P(treated) = P(X >= 1) = 1 - exp(-lambda * R)

    An economic-threshold sprayer treats it only if enough weeds are seen for
    the expected yield loss to exceed the herbicide cost:

        P(treated) = P(X >= m),  m from _et_threshold_count()

    CONSEQUENCES, which are the point of the change:
      * a cell holding few weeds is usually missed, and is missed again the next
        year, so its seedbank compounds instead of being reset;
      * the probability rises steeply with density, so an escaping patch is
        eventually found and treated - escapes grow but are self-limiting;
      * because lambda scales with CELL AREA, the same weed density is harder to
        detect in a 25 x 25 cm cell (0.0625 m2) than in a 1 x 1 m cell. Finer
        spraying resolutions therefore leak more weeds, independently of the
        Year-0 effect already in the model.
    """
    lam_tot = lam_d + lam_m
    if lam_tot <= 0.0:
        return 0.0
    mu = lam_tot * recall              # mean number DETECTED in this cell

    if scenario in ('SCENARIO_2_SPOT', 'SCENARIO_3_SECURE'):
        return 1.0 - math.exp(-mu) if mu < 700.0 else 1.0

    # ET and DynET: need at least m detected weeds
    m = _et_threshold_count(lam_d, lam_m, area, su, scenario)
    if m is None:
        return 0.0
    return _poisson_sf(m - 1, mu)


def _merge_states(states):
    """Collapse trajectories whose seedbanks have converged.

    Without this the state count doubles every year. Keys are the seedbank
    values rounded to MERGE_DIGITS significant figures; merged states carry the
    summed probability and the probability-weighted mean seedbank, so no
    probability mass is lost.
    """
    def sig(x):
        if x <= 0:
            return 0.0
        d = MERGE_DIGITS - int(math.floor(math.log10(x))) - 1
        return round(x, d)

    merged = {}
    for bd, bm, pr in states:
        key = (sig(bd), sig(bm))
        if key in merged:
            obd, obm, opr = merged[key]
            tot = opr + pr
            merged[key] = ((obd * opr + bd * pr) / tot,
                           (obm * opr + bm * pr) / tot, tot)
        else:
            merged[key] = (bd, bm, pr)

    out = list(merged.values())
    if len(out) > MAX_STATES:
        out.sort(key=lambda t: -t[2])
        keep, drop = out[:MAX_STATES], out[MAX_STATES:]
        lost = sum(t[2] for t in drop)
        if lost > 0:                     # redistribute so probabilities sum to 1
            scale = 1.0 + lost / sum(t[2] for t in keep)
            keep = [(bd, bm, pr * scale) for bd, bm, pr in keep]
        out = keep
    return out


def et_decision_scaleup(recall, scenario='SCENARIO_1_ECON'):
    """Factor converting a detected count into an estimated true count.

    Returns 1.0 - i.e. no correction, the threshold acts on what the sensor
    reports - for the economic-threshold strategy when ET_INFORMATION is
    'detected'. The forward-looking strategy (DynET) always applies the
    correction, since it is defined as the informed rule.
    """
    if recall <= 0 or ET_FUTURE_RULE != 'consistent':
        return 1.0
    if scenario == 'SCENARIO_1_ECON' and ET_INFORMATION == 'detected':
        return 1.0
    return 1.0 + ((1.0 / recall) - 1.0) * PATCH_FRACTION


# =====================================================================
# CORE FUNCTION: two-class 10-year forward simulation
#
# =====================================================================
def calculate_10yr_liability_v2(
    surv_dens_d,
    surv_dens_m,
    p_d, p_m,
    yr0_penalty,
    pre_d=None,
    pre_m=None,
    scenario='SCENARIO_1_ECON',
    is_broadcast=False,
    recall=0.41,
    stunt_future=None,
    future_kill=0.95,
    dd_base_kill=None,
    cell_area_m2=1.0          # needed to convert density -> count
):
    """Combined two-class 10-year discounted future weed liability (EUR ha-1).

    Returns an NPV. Divide by ANNUITY_FACTOR for the equivalent annual value
    used in the annualised partial budget (Eq. 6).
    """

    if stunt_future is None:
        stunt_future = STUNTED_SEED_PENALTY
    if dd_base_kill is None:
        dd_base_kill = DD_BASE_KILL

    # -- Year 0 seedbank initialisation (Eq. 8, uses PRE-herbicide density) --
    ref_d = pre_d if (pre_d is not None and pre_d > 0) else surv_dens_d
    if ref_d > 0:
        native_bank_d = ref_d / p_d['emergence']
        yr0_seeds_d = ((p_d['max_seeds'] * surv_dens_d) /
                       (1.0 + p_d['intra_comp'] * surv_dens_d)
                       if surv_dens_d > 0 else 0.0)
        bank_d = native_bank_d * p_d['bank_survival'] + \
            yr0_seeds_d * yr0_penalty * p_d['crop_red']
    else:
        bank_d = 0.0

    ref_m = pre_m if (pre_m is not None and pre_m > 0) else surv_dens_m
    if ref_m > 0:
        native_bank_m = ref_m / p_m['emergence']
        yr0_seeds_m = ((p_m['max_seeds'] * surv_dens_m) /
                       (1.0 + p_m['intra_comp'] * surv_dens_m)
                       if surv_dens_m > 0 else 0.0)
        bank_m = native_bank_m * p_m['bank_survival'] + \
            yr0_seeds_m * yr0_penalty * p_m['crop_red']
    else:
        bank_m = 0.0

    # ---------------------------------------------------------------------
    # YEARS 1-10
    #
    # Under FUTURE_DETECTION = 'expected' this is the original single
    # trajectory. Under 'probabilistic' the cell branches each year into a
    # treated and an untreated continuation, weighted by the probability that
    # the sprayer actually sees the weeds present (see _p_treated). The
    # reported liability is the probability-weighted expectation over all
    # surviving trajectories.
    # ---------------------------------------------------------------------
    su = et_decision_scaleup(recall, scenario)
    probabilistic = (FUTURE_DETECTION == 'probabilistic'
                     and not is_broadcast)

    # each state: (seedbank_dicot, seedbank_monocot, probability)
    states = [(bank_d, bank_m, 1.0)]
    total_npv_cost = 0.0

    for year in range(1, SIM_YEARS + 1):
        discount = 1.0 / ((1.0 + DISCOUNT_RATE) ** year)
        year_cost = 0.0
        next_states = []

        for bd, bm, prob in states:
            emerged_d = bd * p_d['emergence']
            emerged_m = bm * p_m['emergence']

            # -- probability that the cell is treated this year --------------
            if is_broadcast or scenario == 'SCENARIO_3_SECURE':
                # both always pass over the whole field; for DD the BRANCH is
                # full dose vs background dose, not treated vs untreated
                p_treat = 1.0
            elif not probabilistic:
                # 'expected': detected = emerged * recall, above zero whenever
                # any weed emerged, so the decision is deterministic
                det_d = emerged_d * recall * su
                det_m = emerged_m * recall * su
                TCL_det = weed_height_modifier * (
                    det_d * dicot_competition_index +
                    det_m * monocot_competition_index)
                YL_det = multispecies_yield_loss_pct(
                    TCL_det, yield_loss_initial_slope, yield_loss_asymptote_pct)
                loss_ha = (YL_det / 100.0) * weed_free_yield_kg_per_m2 * \
                    10000 * crop_price_eur_per_kg
                if scenario == 'SCENARIO_2_SPOT':
                    fires = (emerged_d + emerged_m) * recall > 0
                elif scenario == 'SCENARIO_4_DYNET':
                    fires = (loss_ha + loss_ha / (1.0 + DISCOUNT_RATE)
                             ) > MAX_HERB_COST_YR
                else:
                    fires = loss_ha > MAX_HERB_COST_YR
                p_treat = 1.0 if fires else 0.0
            else:
                p_treat = _p_treated(emerged_d * cell_area_m2,
                                     emerged_m * cell_area_m2,
                                     cell_area_m2, recall, scenario, su)

            # -- evaluate both continuations ---------------------------------
            for treated, branch_p in ((True, p_treat), (False, 1.0 - p_treat)):
                if branch_p <= 0.0:
                    continue
                w = prob * branch_p
                if w < PRUNE_PROB:
                    continue

                if scenario == 'SCENARIO_3_SECURE' and not is_broadcast:
                    # DD: 'treated' means the cell was SEEN and got the full
                    # dose; otherwise it still receives the background dose.
                    if probabilistic:
                        p_full = _p_treated(emerged_d * cell_area_m2,
                                            emerged_m * cell_area_m2,
                                            cell_area_m2, recall,
                                            'SCENARIO_2_SPOT', su)
                    else:
                        p_full = 1.0 if (emerged_d + emerged_m) * \
                            recall > 0 else 0.0
                    # branch over dose instead of over treatment
                    for full_dose, dose_p in ((True, p_full), (False, 1.0 - p_full)):
                        wd = prob * dose_p
                        if dose_p <= 0.0 or wd < PRUNE_PROB:
                            continue
                        if full_dose:
                            kill_y = future_kill
                            chem = MAX_HERB_COST_YR
                        else:
                            kill_y = dd_base_kill
                            chem = (MAX_HERB_COST_YR if DD_FUTURE_COST_MODE == 'full_rate'
                                    else MAX_HERB_COST_YR * DD_BASE_DOSE)
                        nbd, nbm, cost = _advance_one_year(
                            bd, bm, emerged_d, emerged_m, True, kill_y,
                            stunt_future, chem, p_d, p_m)
                        year_cost += wd * cost
                        next_states.append((nbd, nbm, wd))
                    break      # dose branching replaces treatment branching

                kill_y = future_kill
                chem = MAX_HERB_COST_YR if treated else 0.0
                nbd, nbm, cost = _advance_one_year(
                    bd, bm, emerged_d, emerged_m, treated, kill_y,
                    stunt_future, chem, p_d, p_m)
                year_cost += w * cost
                next_states.append((nbd, nbm, w))

        total_npv_cost += year_cost * discount
        states = _merge_states(next_states) if len(
            next_states) > 1 else next_states
        _BRANCH_STATS['max_states'] = max(
            _BRANCH_STATS['max_states'], len(states))

    return total_npv_cost


def multispecies_yield_loss_pct(TCL, slope, A):
    if TCL <= 0:
        return 0.0
    linear_impact = slope * TCL
    return min(A, linear_impact / (1.0 + (linear_impact / A)))


def year0_revenue(S_d, S_m, area, sprayed):
    """Year-0 revenue after weed competition.

    Kept as a separate function so revenue can also be recomputed inside the
    recall and background-efficacy sweeps."""
    TCL_surv = weed_height_modifier * ((S_d * dicot_competition_index +
                                        S_m * monocot_competition_index) / area)
    if sprayed:
        TCL_surv *= (1.0 / 3.0)
    YL = multispecies_yield_loss_pct(TCL_surv, yield_loss_initial_slope,
                                     yield_loss_asymptote_pct)
    return (weed_free_yield_kg_per_m2 * (1.0 - YL / 100.0)) * area * crop_price_eur_per_kg


def ensure_field_exists(layer, name, qtype):
    if layer.fields().indexFromName(name) == -1:
        layer.dataProvider().addAttributes([QgsField(name, qtype)])
        layer.updateFields()


# =====================================================================
# CELL-LEVEL MEMOISATION
#
# Everything below is a deterministic function of (P_d, P_m, area) plus the
# module-level configuration. Distinct states are few (measured: ~28 across
# 160,000 cells at 25 x 25 cm), so caching here replaces ~45 s of computation
# per hectare with a single dict lookup per cell.
#
# Field-level constants (drip rates etc.) are set once before the loop and are
# constant within a run, so they do not need to be part of the key.
# =====================================================================
_CELL_CACHE = {}
_CACHE_STATS = {'hit': 0, 'miss': 0}


def compute_cell(P_d, P_m, area):
    if USE_CELL_CACHE:
        key = (P_d, P_m, round(area, 6))
        cached = _CELL_CACHE.get(key)
        if cached is not None:
            _CACHE_STATS['hit'] += 1
            return cached
    _CACHE_STATS['miss'] += 1
    res = _compute_cell_uncached(P_d, P_m, area)
    if USE_CELL_CACHE:
        _CELL_CACHE[key] = res
    return res


def _compute_cell_uncached(P_d, P_m, area):
    """All per-cell economics. Returns a dict of contributions to accumulate."""
    area_ha = area / 10000.0
    out = {}

    # -- True count estimation, baseline 70/30 split (Eq. 7) ----------------
    T_d = (P_d + P_d * ((1.0 / recall_dicot) - 1.0)
           * 0.70) if P_d > 0 else isolated_drip_d
    T_m = (P_m + P_m * ((1.0 / recall_monocot) - 1.0)
           * 0.70) if P_m > 0 else isolated_drip_m

    bio_dens_d = T_d / area
    bio_dens_m = T_m / area

    # -- Year 0 spray decision ----------------------------------------------
    # The DECISION density depends on what information the rule is
    # assumed to act on. Under ET_INFORMATION = 'detected' the economic
    # threshold sees only the raw detected counts, as a commercial system would;
    # under 'true' it sees the recall-corrected estimate. The BIOLOGY always
    # uses the estimated true counts (bio_dens_*): the weeds that are present
    # cause yield loss and set seed whether or not the sensor found them.
    if SCENARIO == 'SCENARIO_1_ECON' and ET_INFORMATION == 'detected':
        dec_dens_d, dec_dens_m = P_d / area, P_m / area
    else:
        dec_dens_d, dec_dens_m = bio_dens_d, bio_dens_m

    map_TCL = weed_height_modifier * (dec_dens_d * dicot_competition_index +
                                      dec_dens_m * monocot_competition_index)
    YL_decision = multispecies_yield_loss_pct(map_TCL, yield_loss_initial_slope,
                                              yield_loss_asymptote_pct)
    loss_eur_per_m2 = (YL_decision / 100.0) * \
        weed_free_yield_kg_per_m2 * crop_price_eur_per_kg

    spray = False
    dose = 0.0
    kill = 0.0

    if SCENARIO == 'SCENARIO_1_ECON':
        if (loss_eur_per_m2 > herbicide_cost_fullrate_eur_per_m2) and ((P_d + P_m) > 0):
            spray = True
            dose = 1.0

    elif SCENARIO == 'SCENARIO_2_SPOT':
        if (P_d + P_m) > 0:
            spray = True
            dose = 1.0

    elif SCENARIO == 'SCENARIO_3_SECURE':
        spray = True
        dose = 1.0 if (P_d + P_m) > 0 else DD_BASE_DOSE

    elif SCENARIO == 'SCENARIO_4_DYNET':
        # FORWARD-LOOKING ECONOMIC THRESHOLD
        #
        # Compare full economics of both branches instead of Year-0 yield loss
        # alone. This internalises the shadow cost of letting the cell set seed
        # without needing it as a separate input:
        #
        #   spray if (rev_spray - herb - liab_spray) > (rev_nospray - liab_nospray)
        #
        # As with ET, a cell can only be sprayed if at least one weed was
        # actually DETECTED — you cannot target a cell you believe is empty.
        if (P_d + P_m) > 0:
            S_d_s = T_d * (1.0 - full_kill_rate)
            S_m_s = T_m * (1.0 - full_kill_rate)
            rev_s = year0_revenue(S_d_s, S_m_s, area, True)
            liab_s = calculate_10yr_liability_v2(
                S_d_s / area, S_m_s / area, p_dicot, p_mono, STUNTED_SEED_PENALTY,
                pre_d=bio_dens_d, pre_m=bio_dens_m,
                scenario=SCENARIO, recall=recall_dicot, cell_area_m2=area) * area_ha / ANNUITY_FACTOR

            rev_n = year0_revenue(T_d, T_m, area, False)
            liab_n = calculate_10yr_liability_v2(
                T_d / area, T_m / area, p_dicot, p_mono, 1.0,
                pre_d=bio_dens_d, pre_m=bio_dens_m,
                scenario=SCENARIO, recall=recall_dicot, cell_area_m2=area) * area_ha / ANNUITY_FACTOR

            herb_if_spray = herbicide_cost_fullrate_eur_per_m2 * area
            if (rev_s - herb_if_spray - liab_s) > (rev_n - liab_n):
                spray = True
                dose = 1.0

    if spray:
        kill = full_kill_rate if dose == 1.0 else DD_BASE_KILL
        real_herb = (herbicide_cost_fullrate_eur_per_m2 * area) * dose
        yr0_penalty = STUNTED_SEED_PENALTY
    else:
        kill = 0.0
        dose = 0.0
        real_herb = 0.0
        yr0_penalty = 1.0

    S_d = T_d * (1.0 - kill)
    S_m = T_m * (1.0 - kill)
    dens_d = S_d / area
    dens_m = S_m / area

    # -- FUTURE LIABILITY, SSWM ---------------------------------------------
    shad_sswm_base = calculate_10yr_liability_v2(
        dens_d, dens_m, p_dicot, p_mono, yr0_penalty,
        pre_d=bio_dens_d, pre_m=bio_dens_m,
        scenario=SCENARIO, recall=recall_dicot, cell_area_m2=area) * area_ha

    shad_sswm_low = calculate_10yr_liability_v2(
        dens_d, dens_m, p_dicot_low, p_mono_low, yr0_penalty,
        pre_d=bio_dens_d, pre_m=bio_dens_m,
        scenario=SCENARIO, recall=recall_dicot, cell_area_m2=area) * area_ha

    shad_sswm_high = calculate_10yr_liability_v2(
        dens_d, dens_m, p_dicot_high, p_mono_high, yr0_penalty,
        pre_d=bio_dens_d, pre_m=bio_dens_m,
        scenario=SCENARIO, recall=recall_dicot, cell_area_m2=area) * area_ha

    # -- MARGINAL SHADOW COST ------------------------------------------------
    npv_plus_d = calculate_10yr_liability_v2(
        (T_d + 1.0) * (1.0 - kill) / area, dens_m,
        p_dicot, p_mono, yr0_penalty,
        pre_d=bio_dens_d + (1.0 / area), pre_m=bio_dens_m,
        scenario=SCENARIO, recall=recall_dicot, cell_area_m2=area) * area_ha
    marg_sc_d = npv_plus_d - shad_sswm_base

    npv_plus_m = calculate_10yr_liability_v2(
        dens_d, (T_m + 1.0) * (1.0 - kill) / area,
        p_dicot, p_mono, yr0_penalty,
        pre_d=bio_dens_d, pre_m=bio_dens_m + (1.0 / area),
        scenario=SCENARIO, recall=recall_dicot, cell_area_m2=area) * area_ha
    marg_sc_m = npv_plus_m - shad_sswm_base

    if ANNUALISE_SHADOW_COST:
        marg_sc_d /= ANNUITY_FACTOR
        marg_sc_m /= ANNUITY_FACTOR

    # -- SPLIT SENSITIVITY, SSWM --------------------------------------------
    split_sswm = {}
    for split_name, (patch_frac, drip_d_var, drip_m_var) in drip_rates.items():
        T_d_v = (P_d + P_d * ((1.0 / recall_dicot) - 1.0)
                 * patch_frac) if P_d > 0 else drip_d_var
        T_m_v = (P_m + P_m * ((1.0 / recall_monocot) - 1.0)
                 * patch_frac) if P_m > 0 else drip_m_var
        split_sswm[split_name] = calculate_10yr_liability_v2(
            T_d_v * (1.0 - kill) / area, T_m_v * (1.0 - kill) / area,
            p_dicot, p_mono, yr0_penalty,
            pre_d=T_d_v / area, pre_m=T_m_v / area,
            scenario=SCENARIO, recall=recall_dicot, cell_area_m2=area) * area_ha

    # -- KILL SENSITIVITY, SSWM ---------------------------------------------
    kill_sswm = {}
    for k_name, k_params in KILL_SCENARIOS.items():
        # Year-0 kill must respect the applied DOSE. Under DD, a
        # background-dose cell kills at DD_BASE_KILL, not at the full-rate
        # efficacy being varied here; the background efficacy has its own
        # dedicated sensitivity (DD_KILL_SWEEP -> Fig. 6). Without this guard
        # the k_base column would no longer reproduce the DD baseline once
        # DD_BASE_KILL != full_kill_rate. ET/PB/DynET are unaffected: they
        # always spray at dose == 1.0. SELF-CHECK: for every strategy,
        # Kill_Base_Liab must equal the baseline liability of the same run.
        kill_yr0 = ((k_params['kill'] if dose == 1.0 else DD_BASE_KILL)
                    if spray else 0.0)
        stunt_yr0 = k_params['stunt'] if spray else 1.0
        kill_sswm[k_name] = calculate_10yr_liability_v2(
            T_d * (1.0 - kill_yr0) / area, T_m * (1.0 - kill_yr0) / area,
            p_dicot, p_mono, stunt_yr0,
            pre_d=bio_dens_d, pre_m=bio_dens_m,
            scenario=SCENARIO, recall=recall_dicot,
            stunt_future=k_params['stunt'],
            future_kill=k_params['kill'], cell_area_m2=area) * area_ha

    # -- RECALL SENSITIVITY (original +/-20%), SSWM -------------------------
    recall_sswm = {}
    for rname, (rv, rdrip_d, rdrip_m) in recall_precomp.items():
        T_d_r = (P_d + P_d * ((1.0 / rv) - 1.0) * 0.70) if P_d > 0 else rdrip_d
        T_m_r = (P_m + P_m * ((1.0 / rv) - 1.0) * 0.70) if P_m > 0 else rdrip_m
        recall_sswm[rname] = calculate_10yr_liability_v2(
            T_d_r * (1.0 - kill) / area, T_m_r * (1.0 - kill) / area,
            p_dicot, p_mono, yr0_penalty,
            pre_d=T_d_r / area, pre_m=T_m_r / area,
            scenario=SCENARIO, recall=rv, cell_area_m2=area) * area_ha

    # -- RECALL SWEEP: liability AND revenue ------------------------
    sweep_l_sswm = {}
    sweep_r_sswm = {}
    sweep_l_bc = {}
    sweep_r_bc = {}
    if ENABLE_RECALL_SWEEP:
        for rv in RECALL_SWEEP:
            rdrip_d, rdrip_m = sweep_precomp[rv]
            T_d_r = (P_d + P_d * ((1.0 / rv) - 1.0)
                     * 0.70) if P_d > 0 else rdrip_d
            T_m_r = (P_m + P_m * ((1.0 / rv) - 1.0)
                     * 0.70) if P_m > 0 else rdrip_m

            S_d_r = T_d_r * (1.0 - kill)
            S_m_r = T_m_r * (1.0 - kill)
            sweep_r_sswm[rv] = year0_revenue(S_d_r, S_m_r, area, spray)
            sweep_l_sswm[rv] = calculate_10yr_liability_v2(
                S_d_r / area, S_m_r / area, p_dicot, p_mono, yr0_penalty,
                pre_d=T_d_r / area, pre_m=T_m_r / area,
                scenario=SCENARIO, recall=rv, cell_area_m2=area) * area_ha

            S_d_b = T_d_r * (1.0 - full_kill_rate)
            S_m_b = T_m_r * (1.0 - full_kill_rate)
            sweep_r_bc[rv] = year0_revenue(S_d_b, S_m_b, area, True)
            sweep_l_bc[rv] = calculate_10yr_liability_v2(
                S_d_b / area, S_m_b / area, p_dicot, p_mono, STUNTED_SEED_PENALTY,
                pre_d=T_d_r / area, pre_m=T_m_r / area,
                is_broadcast=True, recall=rv, cell_area_m2=area) * area_ha

    # -- DUAL-DOSE BACKGROUND KILL SWEEP ----------------------------
    ddk_l = {}
    ddk_r = {}
    ddk_h = {}
    if ENABLE_DD_KILL_SWEEP and SCENARIO == 'SCENARIO_3_SECURE':
        for ddk in DD_KILL_SWEEP:
            kill_dd = full_kill_rate if dose == 1.0 else ddk
            S_d_dd = T_d * (1.0 - kill_dd)
            S_m_dd = T_m * (1.0 - kill_dd)
            ddk_r[ddk] = year0_revenue(S_d_dd, S_m_dd, area, True)
            ddk_h[ddk] = (herbicide_cost_fullrate_eur_per_m2 * area) * dose
            ddk_l[ddk] = calculate_10yr_liability_v2(
                S_d_dd / area, S_m_dd / area, p_dicot, p_mono, STUNTED_SEED_PENALTY,
                pre_d=bio_dens_d, pre_m=bio_dens_m,
                scenario=SCENARIO, recall=recall_dicot,
                dd_base_kill=ddk, cell_area_m2=area) * area_ha

    # -- CURRENT-YEAR REVENUE (SSWM) ----------------------------------------
    rev_sswm = year0_revenue(S_d, S_m, area, spray)

    # -- PATH B: BROADCAST ---------------------------------------------------
    bc_cost = herbicide_cost_fullrate_eur_per_m2 * area
    dens_d_bc = T_d * (1.0 - full_kill_rate) / area
    dens_m_bc = T_m * (1.0 - full_kill_rate) / area

    shad_bc_base = calculate_10yr_liability_v2(
        dens_d_bc, dens_m_bc, p_dicot, p_mono, STUNTED_SEED_PENALTY,
        pre_d=bio_dens_d, pre_m=bio_dens_m,
        is_broadcast=True, recall=recall_dicot, cell_area_m2=area) * area_ha
    shad_bc_low = calculate_10yr_liability_v2(
        dens_d_bc, dens_m_bc, p_dicot_low, p_mono_low, STUNTED_SEED_PENALTY,
        pre_d=bio_dens_d, pre_m=bio_dens_m,
        is_broadcast=True, recall=recall_dicot, cell_area_m2=area) * area_ha
    shad_bc_high = calculate_10yr_liability_v2(
        dens_d_bc, dens_m_bc, p_dicot_high, p_mono_high, STUNTED_SEED_PENALTY,
        pre_d=bio_dens_d, pre_m=bio_dens_m,
        is_broadcast=True, recall=recall_dicot, cell_area_m2=area) * area_ha

    split_bc = {}
    for split_name, (patch_frac, drip_d_var, drip_m_var) in drip_rates.items():
        T_d_v = (P_d + P_d * ((1.0 / recall_dicot) - 1.0)
                 * patch_frac) if P_d > 0 else drip_d_var
        T_m_v = (P_m + P_m * ((1.0 / recall_monocot) - 1.0)
                 * patch_frac) if P_m > 0 else drip_m_var
        split_bc[split_name] = calculate_10yr_liability_v2(
            T_d_v * (1.0 - full_kill_rate) / area,
            T_m_v * (1.0 - full_kill_rate) / area,
            p_dicot, p_mono, STUNTED_SEED_PENALTY,
            pre_d=T_d_v / area, pre_m=T_m_v / area,
            is_broadcast=True, recall=recall_dicot, cell_area_m2=area) * area_ha

    kill_bc = {}
    for k_name, k_params in KILL_SCENARIOS.items():
        kill_bc[k_name] = calculate_10yr_liability_v2(
            T_d * (1.0 - k_params['kill']) / area,
            T_m * (1.0 - k_params['kill']) / area,
            p_dicot, p_mono, k_params['stunt'],
            pre_d=bio_dens_d, pre_m=bio_dens_m,
            is_broadcast=True, recall=recall_dicot,
            stunt_future=k_params['stunt'],
            future_kill=k_params['kill'], cell_area_m2=area) * area_ha

    recall_bc = {}
    for rname, (rv, rdrip_d, rdrip_m) in recall_precomp.items():
        T_d_r = (P_d + P_d * ((1.0 / rv) - 1.0) * 0.70) if P_d > 0 else rdrip_d
        T_m_r = (P_m + P_m * ((1.0 / rv) - 1.0) * 0.70) if P_m > 0 else rdrip_m
        recall_bc[rname] = calculate_10yr_liability_v2(
            T_d_r * (1.0 - full_kill_rate) / area,
            T_m_r * (1.0 - full_kill_rate) / area,
            p_dicot, p_mono, STUNTED_SEED_PENALTY,
            pre_d=T_d_r / area, pre_m=T_m_r / area,
            is_broadcast=True, recall=rv, cell_area_m2=area) * area_ha

    rev_bc = year0_revenue(T_d * (1.0 - full_kill_rate),
                           T_m * (1.0 - full_kill_rate), area, True)

    if ENABLE_DD_KILL_SWEEP and SCENARIO == 'SCENARIO_3_SECURE':
        ddk_bc = {ddk: shad_bc_base for ddk in DD_KILL_SWEEP}
    else:
        ddk_bc = {}

    out.update(
        T_d=T_d, T_m=T_m, bio_dens_d=bio_dens_d, bio_dens_m=bio_dens_m,
        S_d=S_d, S_m=S_m, dens_d=dens_d, dens_m=dens_m,
        spray=spray, dose=dose, real_herb=real_herb,
        rev_sswm=rev_sswm, rev_bc=rev_bc, bc_cost=bc_cost,
        shad_sswm_base=shad_sswm_base, shad_sswm_low=shad_sswm_low,
        shad_sswm_high=shad_sswm_high,
        shad_bc_base=shad_bc_base, shad_bc_low=shad_bc_low,
        shad_bc_high=shad_bc_high,
        marg_sc_d=marg_sc_d, marg_sc_m=marg_sc_m,
        split_sswm=split_sswm, split_bc=split_bc,
        kill_sswm=kill_sswm, kill_bc=kill_bc,
        recall_sswm=recall_sswm, recall_bc=recall_bc,
        sweep_l_sswm=sweep_l_sswm, sweep_r_sswm=sweep_r_sswm,
        sweep_l_bc=sweep_l_bc, sweep_r_bc=sweep_r_bc,
        ddk_l=ddk_l, ddk_r=ddk_r, ddk_h=ddk_h, ddk_bc=ddk_bc,
    )
    return out


# ================= MAIN EXECUTION =================
lyr = QgsProject.instance().mapLayersByName(GRID_LAYER_NAME)[0]

out_fields = [
    ('Area_m2',   QVariant.Double), ('P_Dicot',   QVariant.Double),
    ('P_Monocot', QVariant.Double), ('TDicot',    QVariant.Double),
    ('TMonocot',  QVariant.Double), ('TrueDenD',  QVariant.Double),
    ('TrueDenM',  QVariant.Double), ('SDicot',    QVariant.Double),
    ('SMonocot',  QVariant.Double), ('SurvDenD',  QVariant.Double),
    ('SurvDenM',  QVariant.Double), ('Spray_Dec', QVariant.String),
    ('Dose_Rate', QVariant.Double), ('Rev_Eur',   QVariant.Double),
    ('Herb_Cost', QVariant.Double), ('Fut_Liab',  QVariant.Double),
    ('Net_Ret',   QVariant.Double), ('Marg_SC_D', QVariant.Double),
    ('Marg_SC_M', QVariant.Double)
]
for nm, t in out_fields:
    ensure_field_exists(lyr, nm, t)

fi = {f.name(): i for i, f in enumerate(lyr.fields())}
idx_d = lyr.fields().indexFromName(INPUT_DICOT_FIELD)
idx_m = lyr.fields().indexFromName(INPUT_MONOCOT_FIELD)

# -- / single pass for all field-level totals -------------
# A single pass is made over the layer, rather than four separate
# getFeatures() passes. Both removed.
print("Loading features and computing field totals (single pass)...")
features_dict = {}
total_detected_d = 0.0
total_detected_m = 0.0
empty_cells_d = 0
empty_cells_m = 0
empty_cells_both = 0          # cells with NO detection of either class

for f in lyr.getFeatures():
    features_dict[f.id()] = f
    a = f.attributes()
    pd_ = float(a[idx_d] or 0.0)
    pm_ = float(a[idx_m] or 0.0)
    total_detected_d += pd_
    total_detected_m += pm_
    if pd_ == 0:
        empty_cells_d += 1
    if pm_ == 0:
        empty_cells_m += 1
    if pd_ == 0 and pm_ == 0:
        empty_cells_both += 1

total_missed_d = (total_detected_d / recall_dicot -
                  total_detected_d) if recall_dicot > 0 else 0.0
total_missed_m = (total_detected_m / recall_monocot -
                  total_detected_m) if recall_monocot > 0 else 0.0

drip_rates = {}
for split_name, patch_frac in SPLIT_SCENARIOS.items():
    iso_frac = 1.0 - patch_frac
    drip_rates[split_name] = (
        patch_frac,
        (total_missed_d * iso_frac) /
        empty_cells_d if empty_cells_d > 0 else 0.0,
        (total_missed_m * iso_frac) /
        empty_cells_m if empty_cells_m > 0 else 0.0,
    )

isolated_drip_d = drip_rates['split_70_30'][1]
isolated_drip_m = drip_rates['split_70_30'][2]

recall_precomp = {}
for rname, rv in RECALL_SCENARIOS.items():
    r_miss_d = (total_detected_d / rv - total_detected_d) if rv > 0 else 0.0
    r_miss_m = (total_detected_m / rv - total_detected_m) if rv > 0 else 0.0
    recall_precomp[rname] = (
        rv,
        (r_miss_d * 0.30) / empty_cells_d if empty_cells_d > 0 else 0.0,
        (r_miss_m * 0.30) / empty_cells_m if empty_cells_m > 0 else 0.0,
    )

sweep_precomp = {}
for rv in RECALL_SWEEP:
    r_miss_d = (total_detected_d / rv - total_detected_d) if rv > 0 else 0.0
    r_miss_m = (total_detected_m / rv - total_detected_m) if rv > 0 else 0.0
    sweep_precomp[rv] = (
        (r_miss_d * 0.30) / empty_cells_d if empty_cells_d > 0 else 0.0,
        (r_miss_m * 0.30) / empty_cells_m if empty_cells_m > 0 else 0.0,
    )

print(
    f"Baseline drip: {isolated_drip_d:.4f} dicots/cell | {isolated_drip_m:.4f} monocots/cell")
print(f"Annuity factor: {ANNUITY_FACTOR:.4f}")
print(f"Cells: {len(features_dict)} | scenario={SCENARIO} | cache={'ON' if USE_CELL_CACHE else 'OFF'}")
print(f"Future detection: {FUTURE_DETECTION} | ET future rule: {ET_FUTURE_RULE} "
      f"| DD future cost: {DD_FUTURE_COST_MODE} | DD_BASE_KILL: {DD_BASE_KILL}")


def equilibrium_check(start_density=13.4, years=100, verbose=True):
    """CALIBRATION DIAGNOSTIC - printed at the start of every run.

    Projects the weed population forward under CONTINUOUS BROADCAST application
    at the baseline kill rate, which is the management these fields have
    actually received, and reports where it settles. If the model is calibrated,
    that equilibrium should sit close to the observed field densities. If it
    runs away, the ten-year projection is not representing the system and the
    seed-set parameter needs revisiting.

    Report this number in the manuscript: it is the check a referee will run.
    """
    pdx = p_dicot
    bank = start_density / pdx['emergence'] * pdx['bank_survival']
    for _ in range(years):
        em = bank * pdx['emergence']
        surv = em * (1.0 - full_kill_rate)
        seeds = (((pdx['max_seeds'] * surv) / (1.0 + pdx['intra_comp'] * surv))
                 * pdx['crop_red'] * STUNTED_SEED_PENALTY)
        bank = bank * pdx['bank_survival'] + seeds
    eq = bank * pdx['emergence']
    if verbose:
        print(
            f"[CALIBRATION] seed set of survivors (gamma) = {STUNTED_SEED_PENALTY}")
        print(f"[CALIBRATION] equilibrium under continuous broadcast: "
              f"{eq:.1f} plants m-2")
        if eq > 60:
            print("              WARNING: the population runs away even under "
                  "broadcast control. The projection is not representing the "
                  "system; lower STUNTED_SEED_PENALTY.")
        elif eq < 1:
            print("              WARNING: the population is driven to extinction "
                  "under broadcast control. Raise STUNTED_SEED_PENALTY.")
    return eq


equilibrium_check()
print("Starting simulation ...")

# Accumulators
field_area_m2 = 0.0
sprayed_area_m2 = 0.0
vol_herb_sswm = 0.0
vol_herb_bc = 0.0
sum_rev_sswm = 0.0
sum_herb_sswm = 0.0
sum_rev_bc = 0.0
sum_herb_bc = 0.0

sum_shad_sswm_base = 0.0
sum_shad_bc_base = 0.0
sum_shad_sswm_low = 0.0
sum_shad_bc_low = 0.0
sum_shad_sswm_high = 0.0
sum_shad_bc_high = 0.0

sum_split_sswm = {n: 0.0 for n in SPLIT_SCENARIOS}
sum_split_bc = {n: 0.0 for n in SPLIT_SCENARIOS}
sum_kill_sswm = {n: 0.0 for n in KILL_SCENARIOS}
sum_kill_bc = {n: 0.0 for n in KILL_SCENARIOS}
sum_recall_sswm = {n: 0.0 for n in RECALL_SCENARIOS}
sum_recall_bc = {n: 0.0 for n in RECALL_SCENARIOS}

sweep_liab_sswm = {r: 0.0 for r in RECALL_SWEEP}
sweep_liab_bc = {r: 0.0 for r in RECALL_SWEEP}
sweep_rev_sswm = {r: 0.0 for r in RECALL_SWEEP}
sweep_rev_bc = {r: 0.0 for r in RECALL_SWEEP}

ddk_liab_sswm = {k: 0.0 for k in DD_KILL_SWEEP}
ddk_liab_bc = {k: 0.0 for k in DD_KILL_SWEEP}
ddk_rev_sswm = {k: 0.0 for k in DD_KILL_SWEEP}
ddk_herb_sswm = {k: 0.0 for k in DD_KILL_SWEEP}

bulk_updates = {}

# =====================================================================
# MAIN CELL LOOP — now just: compute (cached) then accumulate
# =====================================================================
for feat_id, feat in features_dict.items():
    area = feat.geometry().area()
    if area <= 0:
        continue

    P_d = float(feat.attributes()[idx_d] or 0.0)
    P_m = float(feat.attributes()[idx_m] or 0.0)

    c = compute_cell(P_d, P_m, area)

    field_area_m2 += area
    if c['spray']:
        sprayed_area_m2 += area
    vol_herb_sswm += c['dose'] * area
    vol_herb_bc += area

    sum_rev_sswm += c['rev_sswm']
    sum_herb_sswm += c['real_herb']
    sum_rev_bc += c['rev_bc']
    sum_herb_bc += c['bc_cost']

    sum_shad_sswm_base += c['shad_sswm_base']
    sum_shad_sswm_low += c['shad_sswm_low']
    sum_shad_sswm_high += c['shad_sswm_high']
    sum_shad_bc_base += c['shad_bc_base']
    sum_shad_bc_low += c['shad_bc_low']
    sum_shad_bc_high += c['shad_bc_high']

    for n in SPLIT_SCENARIOS:
        sum_split_sswm[n] += c['split_sswm'][n]
        sum_split_bc[n] += c['split_bc'][n]
    for n in KILL_SCENARIOS:
        sum_kill_sswm[n] += c['kill_sswm'][n]
        sum_kill_bc[n] += c['kill_bc'][n]
    for n in RECALL_SCENARIOS:
        sum_recall_sswm[n] += c['recall_sswm'][n]
        sum_recall_bc[n] += c['recall_bc'][n]

    if ENABLE_RECALL_SWEEP:
        for rv in RECALL_SWEEP:
            sweep_liab_sswm[rv] += c['sweep_l_sswm'][rv]
            sweep_rev_sswm[rv] += c['sweep_r_sswm'][rv]
            sweep_liab_bc[rv] += c['sweep_l_bc'][rv]
            sweep_rev_bc[rv] += c['sweep_r_bc'][rv]

    if ENABLE_DD_KILL_SWEEP and SCENARIO == 'SCENARIO_3_SECURE':
        for k in DD_KILL_SWEEP:
            ddk_liab_sswm[k] += c['ddk_l'][k]
            ddk_rev_sswm[k] += c['ddk_r'][k]
            ddk_herb_sswm[k] += c['ddk_h'][k]
            ddk_liab_bc[k] += c['ddk_bc'][k]

    bulk_updates[feat_id] = {
        fi['Area_m2']:   area,             fi['P_Dicot']:   P_d,
        fi['P_Monocot']: P_m,              fi['TDicot']:    c['T_d'],
        fi['TMonocot']:  c['T_m'],         fi['TrueDenD']:  c['bio_dens_d'],
        fi['TrueDenM']:  c['bio_dens_m'],  fi['SDicot']:    c['S_d'],
        fi['SMonocot']:  c['S_m'],         fi['SurvDenD']:  c['dens_d'],
        fi['SurvDenM']:  c['dens_m'],
        fi['Spray_Dec']: 'Spray' if c['spray'] else 'No',
        fi['Dose_Rate']: c['dose'],        fi['Rev_Eur']:   c['rev_sswm'],
        fi['Herb_Cost']: c['real_herb'],   fi['Fut_Liab']:  c['shad_sswm_base'],
        fi['Net_Ret']:   (c['rev_sswm'] - c['real_herb'] - c['shad_sswm_base']),
        fi['Marg_SC_D']: c['marg_sc_d'],   fi['Marg_SC_M']: c['marg_sc_m']
    }

print("Writing results to map...")
lyr.dataProvider().changeAttributeValues(bulk_updates)
lyr.triggerRepaint()
iface.mapCanvas().refresh()

hits = _CACHE_STATS['hit']
misses = _CACHE_STATS['miss']
tot = hits + misses
print(f"[PERF] cell cache: {hits} hits / {misses} full computations "
      f"({100.0*hits/tot if tot else 0:.2f}% hit rate)")
if FUTURE_DETECTION == 'probabilistic':
    print(f"[PERF] max simultaneous trajectories in any cell-year: "
          f"{_BRANCH_STATS['max_states']} (ceiling {MAX_STATES})")
    if _BRANCH_STATS['max_states'] >= MAX_STATES:
        print("       WARNING: hit the ceiling. Raise MAX_STATES or lower "
              "MERGE_DIGITS if you want a finer expectation.")

# ================= SUMMARY =================
annuity_factor = ANNUITY_FACTOR
field_area_ha = field_area_m2 / 10000.0


def calc_extra_sc(shad_sswm, shad_bc):
    """Annual-equivalent EXTRA liability of SSWM relative to broadcast.
    The NPV difference is converted with the capital recovery factor
    (1/ANNUITY_FACTOR). THIS IS THE STEP MISSING FROM EQ. 11 IN THE PAPER."""
    return ((shad_sswm - shad_bc) / field_area_ha) / annuity_factor


sc_extra_base = calc_extra_sc(sum_shad_sswm_base, sum_shad_bc_base)
sc_extra_low = calc_extra_sc(sum_shad_sswm_low,  sum_shad_bc_low)
sc_extra_high = calc_extra_sc(sum_shad_sswm_high, sum_shad_bc_high)

net_ann_sswm = ((sum_rev_sswm / field_area_ha) - (sum_herb_sswm / field_area_ha)
                - sc_extra_base - sswm_annuity_per_ha)
net_ann_bc = (sum_rev_bc / field_area_ha) - (sum_herb_bc / field_area_ha)
net_sswm_low = ((sum_rev_sswm / field_area_ha) - (sum_herb_sswm / field_area_ha)
                - sc_extra_low - sswm_annuity_per_ha)
net_sswm_high = ((sum_rev_sswm / field_area_ha) - (sum_herb_sswm / field_area_ha)
                 - sc_extra_high - sswm_annuity_per_ha)

sprayed_pct = (sprayed_area_m2 / field_area_m2) * \
    100 if field_area_m2 > 0 else 0
vol_reduct_pct = (1.0 - vol_herb_sswm / vol_herb_bc) * \
    100 if vol_herb_bc > 0 else 0

# -- PARTIAL BUDGET DECOMPOSITION -----------------------------------
# The four components of the difference vs broadcast. All were computable in
# but only d_yield was printed, and none were exported.
d_herb_saved = (sum_herb_bc - sum_herb_sswm) / field_area_ha   # (a) positive
d_yield_yr0 = (sum_rev_sswm - sum_rev_bc) / field_area_ha    # (b) negative
d_liab = -sc_extra_base                                  # (c) negative
d_tech = -sswm_annuity_per_ha                            # (d) negative
d_net_check = d_herb_saved + d_yield_yr0 + d_liab + d_tech    # (e)

# Break-even technology cost — needs NO extra simulation.
tech_breakeven = d_herb_saved + d_yield_yr0 + d_liab

print(f"\n{'='*70}")
print(f"RESULTS: {SCENARIO}  |  Tech: {TECH_SCENARIO}")
print(f"Field: {field_area_ha:.2f} ha  |  Sprayed: {sprayed_pct:.1f}%")
print(f"{'='*70}")
print(f"{'METRIC':<24} | {'SSWM':>10} | {'BROADCAST':>10} | {'DIFF':>9}")
print(f"{'-'*70}")
print(f"{'Revenue (EUR/ha)':<24} | {sum_rev_sswm/field_area_ha:>10.2f} | {sum_rev_bc/field_area_ha:>10.2f} | {d_yield_yr0:>+9.2f}")
print(f"{'Herbicide (EUR/ha)':<24} | {sum_herb_sswm/field_area_ha:>10.2f} | {sum_herb_bc/field_area_ha:>10.2f} | {d_herb_saved:>+9.2f}")
print(f"{'Future Liability':<24} | {sc_extra_base:>10.2f} | {'0.00':>10} | {d_liab:>+9.2f}")
print(f"{'Tech Annuity':<24} | {sswm_annuity_per_ha:>10.2f} | {'0.00':>10} | {d_tech:>+9.2f}")
print(f"{'-'*70}")
print(f"{'NET RETURN (EUR/ha)':<24} | {net_ann_sswm:>10.2f} | {net_ann_bc:>10.2f} | {net_ann_sswm-net_ann_bc:>+9.2f}")
print(f"{'='*70}")
print(f"Herbicide volume reduced by {vol_reduct_pct:.1f}%")

print(f"\n===== [NEW] PARTIAL BUDGET DECOMPOSITION (EUR/ha vs broadcast) =====")
print(f"  (a) Herbicide cost saved      : {d_herb_saved:>+8.2f}")
print(
    f"  (b) Year-0 yield loss         : {d_yield_yr0:>+8.2f}   <-- never reported in the paper")
print(f"  (c) Future weed penalty       : {d_liab:>+8.2f}")
print(f"  (d) Technology cost           : {d_tech:>+8.2f}")
print(f"  {'-'*44}")
print(f"  (e) NET vs broadcast          : {d_net_check:>+8.2f}")
print(f"\n  Break-even technology cost    : {tech_breakeven:>8.2f} EUR/ha")
print(f"  (max annualised tech+info cost at which SSWM matches broadcast)")
print(f"{'='*70}\n")

print("===== SENSITIVITY 1: Weed Population Dynamics (+/-20%) =====")
print(
    f"  Optimistic: Liability={sc_extra_low:>7.2f}  Net={net_sswm_low:>8.2f}")
print(
    f"  Baseline:   Liability={sc_extra_base:>7.2f}  Net={net_ann_sswm:>8.2f}")
print(
    f"  Pessimist:  Liability={sc_extra_high:>7.2f}  Net={net_sswm_high:>8.2f}")
print(f"{'='*70}\n")

sc_split = {n: calc_extra_sc(
    sum_split_sswm[n], sum_split_bc[n]) for n in SPLIT_SCENARIOS}
sc_kill = {n: calc_extra_sc(
    sum_kill_sswm[n],  sum_kill_bc[n]) for n in KILL_SCENARIOS}


def net_from_sc(sc):
    return (sum_rev_sswm/field_area_ha) - (sum_herb_sswm/field_area_ha) - sc - sswm_annuity_per_ha


print("===== SENSITIVITY 2: False Negative Spatial Allocation =====")
for sn in SPLIT_SCENARIOS:
    print(
        f"  {sn:<14}: Liability={sc_split[sn]:>7.2f}  Net={net_from_sc(sc_split[sn]):>8.2f}")
print(f"{'='*70}\n")

print("===== SENSITIVITY 3: Herbicide Efficacy =====")
for kn, kp in KILL_SCENARIOS.items():
    print(f"  {kn:<10} (k={kp['kill']}, stunt={kp['stunt']}): "
          f"Liability={sc_kill[kn]:>7.2f}  Net={net_from_sc(sc_kill[kn]):>8.2f}")
print(f"{'='*70}\n")

sc_recall = {n: calc_extra_sc(
    sum_recall_sswm[n], sum_recall_bc[n]) for n in RECALL_SCENARIOS}

print("===== SENSITIVITY 4: Sensor Recall (+/-20%, LIABILITY CHANNEL ONLY) =====")
print("  This legacy sensitivity varies the liability only")
print("  and holds Year-0 revenue at the baseline recall, so it does NOT")
print("  measure the full effect of detection quality. Use the sweep below.")
for rn, rv in RECALL_SCENARIOS.items():
    print(f"  {rn} (recall={rv}, FN={round((1-rv)*100,1)}%): "
          f"Liability={sc_recall[rn]:>7.2f}  Net={net_from_sc(sc_recall[rn]):>8.2f}")
print(f"{'='*70}\n")

sweep_rows = []
if ENABLE_RECALL_SWEEP:
    print(
        "===== [NEW] RECALL SWEEP: full effect (liability + Year-0 yield) =====")
    print("  Spray map held FIXED at observed detections => LOWER bound on the")
    print("  value of better detection, so the break-even recall shown is an")
    print("  UPPER bound on what is actually required.")
    print(f"  {'recall':>7} | {'FN%':>5} | {'d_herb':>8} | {'d_yield':>8} | {'d_liab':>8} | {'NET':>8}")
    print(f"  {'-'*58}")
    for rv in RECALL_SWEEP:
        d_h = d_herb_saved
        d_y = (sweep_rev_sswm[rv] - sweep_rev_bc[rv]) / field_area_ha
        d_l = -calc_extra_sc(sweep_liab_sswm[rv], sweep_liab_bc[rv])
        net = d_h + d_y + d_l - sswm_annuity_per_ha
        sweep_rows.append((rv, d_h, d_y, d_l, net))
        print(f"  {rv:>7.2f} | {100*(1-rv):>5.1f} | {d_h:>+8.2f} | {d_y:>+8.2f} | "
              f"{d_l:>+8.2f} | {net:>+8.2f}")
    print(f"{'='*70}\n")

ddk_rows = []
if ENABLE_DD_KILL_SWEEP and SCENARIO == 'SCENARIO_3_SECURE':
    print("===== [NEW] DUAL-DOSE BACKGROUND KILL SWEEP =====")
    print("  Tests the assumption that a 50% background dose controls weeds as")
    print("  well as a full dose. Under 'full_rate' the k_bg = 0.95 differential is")
    print("  zero by construction; under 'dose_scaled' it keeps the future")
    print(f"  chemical saving. Mode: {DD_FUTURE_COST_MODE}. SELF-CHECK: the")
    print(
        f"  k_bg = {DD_BASE_KILL:.2f} row must equal the baseline of this run.")
    print(
        f"  {'k_bg':>6} | {'d_herb':>8} | {'d_yield':>8} | {'d_liab':>8} | {'NET':>8}")
    print(f"  {'-'*50}")
    for ddk in DD_KILL_SWEEP:
        d_h = (sum_herb_bc - ddk_herb_sswm[ddk]) / field_area_ha
        d_y = (ddk_rev_sswm[ddk] - sum_rev_bc) / field_area_ha
        d_l = -calc_extra_sc(ddk_liab_sswm[ddk], ddk_liab_bc[ddk])
        net = d_h + d_y + d_l - sswm_annuity_per_ha
        ddk_rows.append((ddk, d_h, d_y, d_l, net))
        print(
            f"  {ddk:>6.2f} | {d_h:>+8.2f} | {d_y:>+8.2f} | {d_l:>+8.2f} | {net:>+8.2f}")
    print(f"{'='*70}\n")

# ================= CSV EXPORT =================
project_name = QgsProject.instance().baseName() or "Untitled_Project"

csv_row = [
    project_name, SCENARIO, round(field_area_ha, 2), round(sprayed_pct, 2),
    round(sum_rev_sswm / field_area_ha,
          2), round(sum_herb_sswm / field_area_ha, 2),
    round(sc_extra_base, 2), round(
        sswm_annuity_per_ha, 2), round(net_ann_sswm, 2),
    round(sum_rev_bc / field_area_ha, 2), round(sum_herb_bc / field_area_ha, 2),
    0.00, 0.00, round(net_ann_bc, 2),
    round(vol_reduct_pct, 2),
    round(sc_extra_low, 2),  round(net_sswm_low, 2),
    round(sc_extra_high, 2), round(net_sswm_high, 2),
    round(sc_split['split_90_10'], 2), round(
        net_from_sc(sc_split['split_90_10']), 2),
    round(sc_split['split_70_30'], 2), round(
        net_from_sc(sc_split['split_70_30']), 2),
    round(sc_split['split_50_50'], 2), round(
        net_from_sc(sc_split['split_50_50']), 2),
    round(sc_kill['k_low'],  2), round(net_from_sc(sc_kill['k_low']),  2),
    round(sc_kill['k_base'], 2), round(net_from_sc(sc_kill['k_base']), 2),
    round(sc_kill['k_high'], 2), round(net_from_sc(sc_kill['k_high']), 2),
    round(sc_recall['recall_low'],  2), round(
        net_from_sc(sc_recall['recall_low']),  2),
    round(sc_recall['recall_base'], 2), round(
        net_from_sc(sc_recall['recall_base']), 2),
    round(sc_recall['recall_high'], 2), round(
        net_from_sc(sc_recall['recall_high']), 2),
    round(d_herb_saved, 2), round(d_yield_yr0, 2), round(d_liab, 2),
    round(d_tech, 2), round(d_net_check, 2), round(tech_breakeven, 2),
    # run configuration, so a results file always records the model
    # settings that produced it. Prevents mixing runs from different modes.
    FUTURE_DETECTION, ET_FUTURE_RULE, DD_FUTURE_COST_MODE,
    DD_BASE_KILL, recall_dicot, round(ANNUITY_FACTOR, 4),
    # descriptive statistics for the supplementary field table.
    # Captured here so you never have to re-run just to describe the input data.
    len(features_dict),
    round(total_detected_d / field_area_ha, 1),
    round(total_detected_m / field_area_ha, 1),
    round(100.0 * (1.0 - empty_cells_both / len(features_dict)), 2),
]

CSV_HEADER = [
    'Project', 'Scenario', 'Area_Ha', 'Sprayed_Pct',
    'SSWM_Rev', 'SSWM_Herb', 'SSWM_Liab_Base', 'SSWM_Tech', 'SSWM_Net_Base',
    'BC_Rev', 'BC_Herb', 'BC_Liab', 'BC_Tech', 'BC_Net', 'Vol_Red_Pct',
    'Weed_Low_Liab', 'Weed_Low_Net', 'Weed_High_Liab', 'Weed_High_Net',
    'Split_90_10_Liab', 'Split_90_10_Net',
    'Split_70_30_Liab', 'Split_70_30_Net',
    'Split_50_50_Liab', 'Split_50_50_Net',
    'Kill_Low_Liab', 'Kill_Low_Net',
    'Kill_Base_Liab', 'Kill_Base_Net',
    'Kill_High_Liab', 'Kill_High_Net',
    'Recall_Low_Liab',  'Recall_Low_Net',
    'Recall_Base_Liab', 'Recall_Base_Net',
    'Recall_High_Liab', 'Recall_High_Net',
    'D_Herb_Saved', 'D_Yield_Yr0', 'D_Liab', 'D_Tech', 'D_Net', 'Tech_Breakeven',
    'Future_Detection', 'ET_Future_Rule', 'DD_Future_Cost', 'DD_Base_Kill',
    'Recall', 'Annuity_Factor',
    'N_Cells', 'Det_Dicot_Per_Ha', 'Det_Mono_Per_Ha', 'Pct_Cells_Infested',
]

file_exists = os.path.isfile(RESULTS_CSV_PATH)
try:
    os.makedirs(os.path.dirname(RESULTS_CSV_PATH), exist_ok=True)
    with open(RESULTS_CSV_PATH, 'a', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        if not file_exists:
            writer.writerow(CSV_HEADER)
        writer.writerow(csv_row)
    print(f"[SUCCESS] Results appended to:\n{RESULTS_CSV_PATH}")
except Exception as e:
    print(f"[ERROR] Could not write to CSV: {e}")


def _append_rows(path, header, rows, prefix):
    exists = os.path.isfile(path)
    try:
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, 'a', newline='', encoding='utf-8') as fh:
            w = csv.writer(fh)
            if not exists:
                w.writerow(header)
            for r in rows:
                w.writerow(list(prefix) + [round(v, 4) for v in r])
        print(f"[SUCCESS] appended {len(rows)} rows to {path}")
    except Exception as exc:
        print(f"[ERROR] could not write {path}: {exc}")


if ENABLE_RECALL_SWEEP and sweep_rows:
    _append_rows(
        RESULTS_CSV_PATH.replace('.csv', '_recall_sweep.csv'),
        ['Project', 'Scenario', 'Area_Ha', 'Recall',
            'D_Herb', 'D_Yield', 'D_Liab', 'Net'],
        sweep_rows, (project_name, SCENARIO, round(field_area_ha, 2)))

if ENABLE_DD_KILL_SWEEP and ddk_rows:
    _append_rows(
        RESULTS_CSV_PATH.replace('.csv', '_dd_kill_sweep.csv'),
        ['Project', 'Scenario', 'Area_Ha', 'K_Background',
            'D_Herb', 'D_Yield', 'D_Liab', 'Net'],
        ddk_rows, (project_name, SCENARIO, round(field_area_ha, 2)))
