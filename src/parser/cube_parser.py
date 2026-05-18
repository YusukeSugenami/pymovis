import numpy as np

def load_cube(inp_file):

    with open(inp_file) as f:
        lines = f.readlines()

    natoms = abs(int(lines[2].split()[0]))
    origin = np.array([float(x) for x in lines[2].split()[1:4]])
    nx, vx = (int(lines[3].split()[0]), np.array([float(x) for x in lines[3].split()[1:4]]))
    ny, vy = (int(lines[4].split()[0]), np.array([float(x) for x in lines[4].split()[1:4]]))
    nz, vz = (int(lines[5].split()[0]), np.array([float(x) for x in lines[5].split()[1:4]]))
    grid_vecs = np.vstack([vx, vy, vz])

    atomnos = []
    coords = []
    for i in range(natoms):
        parts = lines[6 + i].split()
        atomnos.append(int((parts[0])))
        coords.append([float(parts[2]), float(parts[3]), float(parts[4])])
    atomnos = np.array(atomnos, dtype=int)
    coords = np.array(coords, dtype=float) # bohr

    if len(lines[6 + natoms].strip().split(" ")) < 6:
        raw_values = " ".join(lines[7 + natoms :]).replace("D", "E").replace("d", "E")
    else:
        raw_values = " ".join(lines[6 + natoms :]).replace("D", "E").replace("d", "E")
    data = np.fromstring(raw_values, sep=" ")
    expected_size = nx * ny * nz
    if data.size != expected_size:
        raise ValueError(
            f"cube grid size mismatch: expected {expected_size} values for shape ({nx}, {ny}, {nz}), got {data.size}"
        )
    mo_cube = np.array(data).reshape((nx, ny, nz), order="C")
    
    cube_info = {
        'file_type' : 'cube',
        'atomnos' : atomnos,
        'coords' : coords,
        'coords_unit' : 'Bohr',
        'origin' : origin,
        'grid_vecs' : grid_vecs,
        'mo_cube' : mo_cube
    }

    return cube_info
