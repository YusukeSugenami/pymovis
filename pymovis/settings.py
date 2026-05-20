ELEMENT_DATA = {
    "H": {"color": (1.00, 1.00, 1.00), "radius": 0.31},
    "He": {"color": (0.85, 1.00, 1.00), "radius": 0.28},
    "Li": {"color": (0.80, 0.50, 1.00), "radius": 1.28},
    "Be": {"color": (0.76, 1.00, 0.00), "radius": 0.96},
    "B": {"color": (1.00, 0.71, 0.71), "radius": 0.84},
    "C": {"color": (0.57, 0.57, 0.57), "radius": 0.76},
    "N": {"color": (0.05, 0.39, 1.00), "radius": 0.71},
    "O": {"color": (1.00, 0.05, 0.05), "radius": 0.66},
    "F": {"color": (0.56, 0.88, 0.31), "radius": 0.57},
    "Ne": {"color": (0.70, 0.89, 0.96), "radius": 0.58},
    "Na": {"color": (0.67, 0.36, 0.95), "radius": 1.66},
    "Mg": {"color": (0.54, 1.00, 0.00), "radius": 1.41},
    "Al": {"color": (0.75, 0.65, 0.65), "radius": 1.21},
    "Si": {"color": (0.94, 0.78, 0.63), "radius": 1.11},
    "P": {"color": (1.00, 0.50, 0.00), "radius": 1.07},
    "S": {"color": (1.00, 1.00, 0.18), "radius": 1.05},
    "Cl": {"color": (0.12, 0.94, 0.12), "radius": 1.02},
    "Ar": {"color": (0.50, 0.82, 0.89), "radius": 1.06},
    "K": {"color": (0.56, 0.25, 0.83), "radius": 2.03},
    "Ca": {"color": (0.24, 1.00, 0.00), "radius": 1.76},
    "Sc": {"color": (0.90, 0.90, 0.90), "radius": 1.70},
    "Ti": {"color": (0.75, 0.76, 0.78), "radius": 1.60},
    "V": {"color": (0.65, 0.65, 0.67), "radius": 1.53},
    "Cr": {"color": (0.54, 0.60, 0.78), "radius": 1.39},
    "Mn": {"color": (0.61, 0.48, 0.78), "radius": 1.39},
    "Fe": {"color": (0.88, 0.40, 0.20), "radius": 1.32},
    "Co": {"color": (0.94, 0.56, 0.63), "radius": 1.26},
    "Ni": {"color": (0.31, 0.82, 0.31), "radius": 1.24},
    "Cu": {"color": (1.00, 0.49, 0.31), "radius": 1.32},
    "Zn": {"color": (0.49, 0.50, 0.69), "radius": 1.22},
    "Ga": {"color": (0.76, 0.56, 0.56), "radius": 1.22},
    "Ge": {"color": (0.40, 0.56, 0.56), "radius": 1.20},
    "As": {"color": (0.74, 0.50, 0.89), "radius": 1.19},
    "Se": {"color": (1.00, 0.63, 0.00), "radius": 1.20},
    "Br": {"color": (0.65, 0.16, 0.16), "radius": 1.20},
    "Kr": {"color": (0.36, 0.72, 0.82), "radius": 1.16},
    "Rb": {"color": (0.44, 0.18, 0.69), "radius": 2.16},
    "Sr": {"color": (0.00, 1.00, 0.00), "radius": 1.91},
    "Y": {"color": (0.58, 1.00, 1.00), "radius": 1.62},
    "Zr": {"color": (0.58, 0.88, 0.88), "radius": 1.48},
    "Nb": {"color": (0.45, 0.76, 0.79), "radius": 1.37},
    "Mo": {"color": (0.33, 0.71, 0.71), "radius": 1.30},
    "Tc": {"color": (0.23, 0.62, 0.62), "radius": 1.27},
    "Ru": {"color": (0.14, 0.56, 0.56), "radius": 1.25},
    "Rh": {"color": (0.04, 0.49, 0.55), "radius": 1.25},
    "Pd": {"color": (0.00, 0.41, 0.52), "radius": 1.20},
    "Ag": {"color": (0.75, 0.75, 0.75), "radius": 1.28},
    "Cd": {"color": (1.00, 0.85, 0.56), "radius": 1.36},
    "In": {"color": (0.65, 0.46, 0.45), "radius": 1.42},
    "Sn": {"color": (0.40, 0.50, 0.50), "radius": 1.39},
    "Sb": {"color": (0.62, 0.39, 0.71), "radius": 1.39},
    "Te": {"color": (0.83, 0.48, 0.00), "radius": 1.38},
    "I": {"color": (0.58, 0.00, 0.58), "radius": 1.39},
    "Xe": {"color": (0.26, 0.62, 0.69), "radius": 1.40},
    "Cs": {"color": (0.34, 0.09, 0.56), "radius": 2.35},
    "Ba": {"color": (0.00, 0.78, 0.00), "radius": 1.98},
    "La": {"color": (0.44, 0.83, 1.00), "radius": 1.69},
}
DEFAULT_DATA = {"color": (0.0, 0.0, 0.0), "radius": 0.7}

ORBITAL_COLORS = {"pos": "blue", "neg": "red"}
ORBITAL_OPACITY = 1.0

BOND_RADIUS = 0.15
BOND_DEFAULT_COLOR = "white"

PADDING = 5.0

IMAGE_SUFFIXES = ("png", "jpeg", "jpg", "bmp", "tif", "tiff")
DEFAULT_OUT_SUFFIX = "png"

ANTI_ALIASING = "ssaa"  # ssaa / msaa / fxaa
BACKGROUND_COLOR = [255, 255, 255, 255]
IMAGE_QUALITY = 4

GRID_SETTING = "shape"  # size / shape
GRID_SHAPE = [100, 100, 100]
GRID_SIZE = [0.3333, 0.3333, 0.3333]

##############################################
# Parser arguments configuration
##############################################

PARSER_ARGS = [
    # Positional arguments
    {
        "args": ["input_file"],
        "kwargs": {
            "help": "Gaussian fchk or cube file to visualize.",
        },
    },
    {
        "args": ["mo"],
        "kwargs": {
            "nargs": "*",
            "default": [-1],
            "help": (
                "MO index to visualize. "
                "Identify mo by number starts from 0 or using keywords (e.g. HOMO, LUMO+2 etc). "
                "Multiple indices can be set"
            ),
        },
    },
    # Optional arguments
    {
        "args": ["-o", "--out"],
        "kwargs": {
            "nargs": "*",
            "default": None,
            "help": "Output file name. Default: {fchk}_{mo}.png",
        },
    },
    {
        "args": ["-i", "--iso"],
        "kwargs": {
            "default": 0.05,
            "type": float,
            "help": "Iso surface value. Default: 0.05",
        },
    },
    {
        "args": ["-b", "--basis"],
        "kwargs": {
            "default": None,
            "help": (
                "Basis function is retrieved from fchk file by default. "
                "If you want to set basis function manually or fchk file does not contain basis information, "
                "use this option like \"-b 6-31G(d,p)\""
            ),
        },
    },
    {
        "args": ["-t", "--transparent"],
        "kwargs": {
            "default": False,
            "type": bool,
            "help": "Transparent background, True / False. Default: False",
        },
    },
    {
        "args": ["-cp", "--camera_pos"],
        "kwargs": {
            "nargs": 3,
            "type": float,
            "default": None,
            "metavar": ("X", "Y", "Z"),
            "help": "Camera position",
        },
    },
    {
        "args": ["-cf", "--camera_focal"],
        "kwargs": {
            "nargs": 3,
            "type": float,
            "default": None,
            "metavar": ("X", "Y", "Z"),
            "help": "Focal point of camera",
        },
    },
    {
        "args": ["-cu", "--camera_up"],
        "kwargs": {
            "nargs": 3,
            "type": float,
            "default": None,
            "metavar": ("X", "Y", "Z"),
            "help": "View-up vector of camera",
        },
    },
    {
        "args": ["-ca", "--camera_axis"],
        "kwargs": {
            "default": None,
            "choices": ["X", "Y", "Z", "x", "y", "z"],
            "help": (
                "Auto camera axis. "
                "If set, camera position is computed to fit the molecule along the selected axis."
            ),
        },
    },
    {
        "args": ["--camera_focal_atoms"],
        "kwargs": {
            "nargs": "+",
            "type": int,
            "default": None,
            "metavar": "ATOM",
            "help": (
                "Atom indices (0-based) used to determine the camera focal point. "
                "One atom uses that atom position; two or more use the centroid."
            ),
        },
    },
    {
        "args": ["--camera_axis_atoms"],
        "kwargs": {
            "nargs": "+",
            "type": int,
            "default": None,
            "metavar": "ATOM",
            "help": (
                "Atom indices (0-based) used to determine the camera view axis. "
                "Two atoms define a line; three or more define a plane normal."
            ),
        },
    },
    {
        "args": ["--camera_up_atoms"],
        "kwargs": {
            "nargs": "+",
            "type": int,
            "default": None,
            "metavar": "ATOM",
            "help": (
                "Atom indices (0-based) used to determine the camera up axis. "
                "Two atoms define a line; three or more define a plane normal."
            ),
        },
    },
    {
        "args": ["--camera_axis_reference"],
        "kwargs": {
            "nargs": "+",
            "type": float,
            "default": None,
            "metavar": "VALUE",
            "help": (
                "Reference vector to fix the sign of the camera view axis. "
                "Can be specified as: "
                "(1) Two atom indices [idx1 idx2] -> vector = coords[idx2] - coords[idx1], or "
                "(2) Three coordinates [X Y Z] -> direct vector. "
                "Used with --camera_axis_atoms; axis direction is flipped "
                "if its dot product with the reference is negative."
            ),
        },
    },
    {
        "args": ["--camera_up_reference"],
        "kwargs": {
            "nargs": "+",
            "type": float,
            "default": None,
            "metavar": "VALUE",
            "help": (
                "Reference vector to fix the sign of the camera up axis. "
                "Can be specified as: "
                "(1) Two atom indices [idx1 idx2] -> vector = coords[idx2] - coords[idx1], or "
                "(2) Three coordinates [X Y Z] -> direct vector. "
                "Used with --camera_up_atoms; axis direction is flipped "
                "if its dot product with the reference is negative."
            ),
        },
    },
    {
        "args": ["-womo", "--without_mo"],
        "kwargs": {
            "action": "store_true",
            "help": "Only render molecule structure without MO isosurface",
        },
    }
]
