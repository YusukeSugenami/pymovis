import numpy as np


AXIS_VIEW_CONFIG = {
    "x": {
        "view": np.array([1.0, 0.0, 0.0]),
        "right": np.array([0.0, 1.0, 0.0]),
        "up": np.array([0.0, 0.0, 1.0]),
    },
    "y": {
        "view": np.array([0.0, 1.0, 0.0]),
        "right": np.array([1.0, 0.0, 0.0]),
        "up": np.array([0.0, 0.0, 1.0]),
    },
    "z": {
        "view": np.array([0.0, 0.0, 1.0]),
        "right": np.array([1.0, 0.0, 0.0]),
        "up": np.array([0.0, 1.0, 0.0]),
    },
}


def _normalize_vector(vector: np.ndarray, label: str) -> np.ndarray:
    norm = np.linalg.norm(vector)
    if norm == 0.0:
        raise ValueError(f"{label} vector must not be zero")
    return vector / norm


def _coerce_atom_indices(atom_indices) -> np.ndarray:
    if isinstance(atom_indices, (int, np.integer)):
        indices = np.array([int(atom_indices)], dtype=int)
    else:
        indices = np.asarray(atom_indices, dtype=int)
    if indices.ndim != 1 or indices.size == 0:
        raise ValueError("atom indices must be a non-empty 1D sequence")
    return indices


def infer_center_from_atom_indices(coords, atom_indices):
    coords = np.asarray(coords, dtype=float)
    indices = _coerce_atom_indices(atom_indices)
    if np.any(indices < 0) or np.any(indices >= coords.shape[0]):
        raise IndexError("atom index out of range")

    selected = coords[indices]
    if indices.size == 1:
        return selected[0]
    return selected.mean(axis=0)


def infer_axis_from_atom_indices(coords, atom_indices):
    coords = np.asarray(coords, dtype=float)
    indices = _coerce_atom_indices(atom_indices)
    if indices.size < 2:
        raise ValueError("atom indices must contain at least two atoms")
    if np.any(indices < 0) or np.any(indices >= coords.shape[0]):
        raise IndexError("atom index out of range")

    selected = coords[indices]
    centered = selected - selected.mean(axis=0)

    if indices.size == 2:
        direction = selected[1] - selected[0]
        return _normalize_vector(direction, "axis")

    # Best-fit plane normal for three or more atoms.
    _, _, vh = np.linalg.svd(centered, full_matrices=False)
    normal = vh[-1]
    return _normalize_vector(normal, "axis")


def _orthonormal_camera_basis(view_axis, up_axis):
    view_axis = _normalize_vector(np.asarray(view_axis, dtype=float), "view")
    up_axis = _normalize_vector(np.asarray(up_axis, dtype=float), "up")

    right_axis = np.cross(view_axis, up_axis)
    right_axis = _normalize_vector(right_axis, "right")
    up_axis = _normalize_vector(np.cross(right_axis, view_axis), "up")

    return view_axis, right_axis, up_axis


def _fit_camera_from_basis(
    coords,
    view_axis,
    up_axis,
    padding: float = 2.0,
    window_size: tuple[int, int] = (1024, 768),
    vertical_fov_deg: float = 30.0,
    safety_factor: float = 1.05,
    min_distance: float = 1.0,
    center: np.ndarray | list[float] | tuple[float, float, float] | None = None,
):
    coords = np.asarray(coords, dtype=float)
    if coords.ndim != 2 or coords.shape[1] != 3 or coords.shape[0] == 0:
        raise ValueError("coords must be a non-empty array with shape (n_atoms, 3)")
    if padding < 0.0:
        raise ValueError("padding must be non-negative")

    width, height = window_size
    if width <= 0 or height <= 0:
        raise ValueError("window_size values must be positive")
    if vertical_fov_deg <= 0.0 or vertical_fov_deg >= 180.0:
        raise ValueError("vertical_fov_deg must be in the range (0, 180)")
    if safety_factor <= 0.0:
        raise ValueError("safety_factor must be positive")
    if min_distance < 0.0:
        raise ValueError("min_distance must be non-negative")

    centroid = coords.mean(axis=0)
    if center is None:
        camera_focal = centroid
    else:
        camera_focal = np.asarray(center, dtype=float)
        if camera_focal.shape != (3,):
            raise ValueError("center must be a 3-element vector")

    view_axis, right_axis, up_axis = _orthonormal_camera_basis(view_axis, up_axis)
    centered = coords - camera_focal

    half_right = np.max(np.abs(centered @ right_axis)) + padding
    half_up = np.max(np.abs(centered @ up_axis)) + padding
    half_view = np.max(np.abs(centered @ view_axis)) + padding

    aspect = float(width) / float(height)
    v_half = np.deg2rad(vertical_fov_deg) / 2.0
    h_half = np.arctan(np.tan(v_half) * aspect)

    fit_dist = half_view + max(
        half_right / np.tan(h_half),
        half_up / np.tan(v_half),
    )
    distance = max(min_distance, fit_dist * safety_factor)

    camera_pos = camera_focal - view_axis * distance
    camera_up = up_axis
    return camera_pos, camera_focal, camera_up


def infer_axis_camera_fit_from_coords(
    coords,
    axis: str,
    padding: float = 2.0,
    window_size: tuple[int, int] = (1024, 768),
    vertical_fov_deg: float = 30.0,
    safety_factor: float = 1.05,
    min_distance: float = 1.0,
    center: np.ndarray | list[float] | tuple[float, float, float] | None = None,
):
    """
    Determine an axis-aligned camera that fits the whole molecule + padding.

    The camera looks from the negative side of the selected axis to the positive
    direction, with screen orientation fixed as:
        - axis="X": right=+Y, up=+Z
        - axis="Y": right=+X, up=+Z
        - axis="Z": right=+X, up=+Y

    Defaults match the current renderer settings in this project:
        - window_size=(1024, 768)  -> aspect ratio 4:3
        - vertical_fov_deg=30.0

    Parameters
    ----------
    coords : array-like
        Atomic coordinates with shape (n_atoms, 3).
    axis : str
        View axis. One of "X", "Y", "Z" (case-insensitive).
    padding : float, optional
        Extra margin added to half-size in right/up/view directions.
    window_size : tuple[int, int], optional
        Render size used to derive aspect ratio (width, height).
    vertical_fov_deg : float, optional
        Vertical field of view in degrees.
    safety_factor : float, optional
        Additional multiplicative margin for robustness.
    min_distance : float, optional
        Lower bound of camera distance from centroid.
    center : array-like, optional
        Custom screen center (camera focal point). If None, centroid is used.

    Returns
    -------
    tuple[np.ndarray, np.ndarray, np.ndarray]
        (camera_pos, camera_focal, camera_up)
    """

    axis_key = str(axis).strip().lower()
    if axis_key not in AXIS_VIEW_CONFIG:
        raise ValueError("axis must be one of 'X', 'Y', or 'Z'")
    basis = AXIS_VIEW_CONFIG[axis_key]
    return _fit_camera_from_basis(
        coords,
        basis["view"],
        basis["up"],
        padding=padding,
        window_size=window_size,
        vertical_fov_deg=vertical_fov_deg,
        safety_factor=safety_factor,
        min_distance=min_distance,
        center=center,
    )


def infer_atom_axis_camera_fit_from_coords(
    coords,
    camera_axis_atoms,
    camera_up_atoms,
    padding: float = 2.0,
    window_size: tuple[int, int] = (1024, 768),
    vertical_fov_deg: float = 30.0,
    safety_factor: float = 1.05,
    min_distance: float = 1.0,
    center: np.ndarray | list[float] | tuple[float, float, float] | None = None,
):
    """
    Determine camera settings from atom-index-defined axis and up directions.

    For each axis option:
        - two atoms define a line axis using the direction from the first to the second atom
        - three or more atoms define a plane axis using the best-fit plane normal
    """

    view_axis = infer_axis_from_atom_indices(coords, camera_axis_atoms)
    up_axis = infer_axis_from_atom_indices(coords, camera_up_atoms)
    return _fit_camera_from_basis(
        coords,
        view_axis,
        up_axis,
        padding=padding,
        window_size=window_size,
        vertical_fov_deg=vertical_fov_deg,
        safety_factor=safety_factor,
        min_distance=min_distance,
        center=center,
    )


