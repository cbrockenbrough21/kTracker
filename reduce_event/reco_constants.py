# reco_constants.py

# Z positions (geometry landmarks)
Z_TARGET = -300.0
Z_DUMP = 42.0

# Sagitta projection tuning constants
SAGITTA_TARGET_CENTER = 1.85
SAGITTA_DUMP_CENTER = 1.5
SAGITTA_TARGET_WIDTH = 0.25
SAGITTA_DUMP_WIDTH = 0.3

# Track slope limit 
TX_MAX = 0.15  # maximum projected slope in X for hodo-to-chamber mapping
TY_MAX = 0.10  # maximum projected slope in Y for hodo-to-chamber mapping
BUFFER = 2  # number of elements to buffer on each side of matched range

N_CHAMBER_PLANES = 30  # IDs 1–30
N_HODO_PLANES    = 8   # IDs 31–38
N_PROP_PLANES    = 16  # IDs 39–54
