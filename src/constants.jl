# Transverse film thicknesses / channel radius; reference: ASETS-II, r=0.5 mm.
const DEFAULT_FILM_MIN_OVER_R = 0.004 # 2 μm at reference radius
const DEFAULT_FILM_TRANSITION_OVER_R = 0.01 # 5 μm at reference radius
const DEFAULT_FILM_MAX_OVER_R = 0.2 # 100 μm at reference radius
const DEFAULT_NUCLEATED_FILM_OVER_R = 0.04 # 20 μm at reference radius

const INCHES = 2.54e-2
const GRAVITY = 9.8
const SEPERATION_VAR = -1e10

const DEFAULT_CLOSED = true
const DEFAULT_BOIL_SCAN_INTERVAL = 1e-2
const DEFAULT_BOIL_WAITING_TIME = 1.0
const DEFAULT_RN_BOIL = 3e-6
const DEFAULT_INERTIA_F = 1.1
const DEFAULT_D = 1e-3
const DEFAULT_TUBESHAPE = "square"
const DEFAULT_NU = 3.6
const DEFAULT_SLUGNUM = 30
const DEFAULT_DFILM_RELATIVE = 0.04
const DEFAULT_FILM_FRACTION = 0.3
const DEFAULT_G = [0.0,0.0]
const DEFAULT_ETAPLUS = 0.6
const DEFAULT_ETAMINUS = 0.0
const DEFAULT_NUCLEATENUM = 250 # total nucleation sites along the whole flow path
# Axial model-resolution scale [m]. Legacy ASETS-II default; override per case.
const DEFAULT_L_NEWBUBBLE = 6e-3
const DEFAULT_LIQUID_CHARGE_RATIO = 0.46
const DEFAULT_SIGMA_CHARGE = 0.01
# Axial lengths / L_newbubble. Names refer to the phases being merged.
const DEFAULT_LIQUID_MERGE_FRAC = 0.5 # short vapor gap -> merge liquids
const DEFAULT_VAPOR_MERGE_FRAC = 2   # short liquid slug -> merge vapors
const DEFAULT_L_ADJUST_MAX_FACTOR = 2 # boiling mass-correction shrink limit
const DEFAULT_L0_FILM_FRAC = 0.15
const DEFAULT_L0_PURE_VAPOR_FRAC = 0.5
const DEFAULT_BOIL_CLEARANCE_FRAC = 2 # from the site to each liquid endpoint
const DEFAULT_NEWBUBBLE_FILM_FRAC = 0.125 # each end of the inserted bubble
