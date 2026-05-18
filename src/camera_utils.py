import numpy as np

def infer_camera_from_coords(coords):

    coords = np.asarray(coords, dtype=float)
    centroid = coords.mean(axis=0)

    if coords.shape[0] >= 3:
        centered = coords - centroid
        cov_matrix = np.cov(centered, rowvar=False)
        eigenvalues, eigenvectors = np.linalg.eigh(cov_matrix)
        view_axis = eigenvectors[:, np.argmin(eigenvalues)]
        up_axis = eigenvectors[:, np.argmax(eigenvalues)]
    elif coords.shape[0] == 2:
        view_axis = np.array([0.0, 0.0, 1.0])
        up_axis = coords[1] - coords[0]
    else:
        view_axis = np.array([0.0, 0.0, 1.0])
        up_axis = np.array([0.0, 1.0, 0.0])

    view_norm = np.linalg.norm(view_axis)
    if view_norm == 0:
        view_axis = np.array([0.0, 0.0, 1.0])
    else:
        view_axis = view_axis / view_norm

    up_norm = np.linalg.norm(up_axis)
    if up_norm == 0:
        up_axis = np.array([0.0, 1.0, 0.0])
    else:
        up_axis = up_axis / up_norm

    extent = np.ptp(coords, axis=0)
    distance = max(float(np.linalg.norm(extent)) * 2.5, 20.0)
    camera_pos = centroid + view_axis * distance

    return camera_pos, centroid, up_axis
