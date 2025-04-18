# Expose classes from submodules for easier import

# Classical Algorithms
try:
    from .classical.dc import DC
except ImportError:
    print("Warning: Could not import DC algorithm.")
    DC = None
try:
    from .classical.bc import BC
except ImportError:
    print("Warning: Could not import BC algorithm.")
    BC = None
try:
    from .classical.cc import CC
except ImportError:
    print("Warning: Could not import CC algorithm.")
    CC = None
try:
    from .classical.ec import EC
except ImportError:
    print("Warning: Could not import EC algorithm.")
    EC = None
try:
    from .classical.ic import IC
except ImportError:
    print("Warning: Could not import IC algorithm.")
    IC = None
try:
    from .classical.nc import NC
except ImportError:
    print("Warning: Could not import NC algorithm.")
    NC = None
try:
    from .classical.sc import SC
except ImportError:
    print("Warning: Could not import SC algorithm.")
    SC = None

# Modern Algorithms
try:
    from .modern.jdc import JDC
except ImportError:
    print("Warning: Could not import JDC algorithm.")
    JDC = None
try:
    from .modern.teo import TEO
except ImportError:
    print("Warning: Could not import TEO algorithm.")
    TEO = None
try:
    from .modern.tgso import TGSO
except ImportError:
    print("Warning: Could not import TGSO algorithm.")
    TGSO = None


# Define __all__ based on successful imports
__all__ = [
    name
    for name, obj in globals().items()
    if obj is not None
    and name
    in [
        # Classical
        "DC",
        "BC",
        "CC",
        "EC",
        "IC",
        "NC",
        "SC",
        # Modern
        "JDC",
        "TEO",
        "TGSO",
    ]
]
