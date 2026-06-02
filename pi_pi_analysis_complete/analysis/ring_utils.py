import numpy as np

def find_all_rings(aromatic_selection):
    ring_atoms = []
    for res in aromatic_selection.residues:
        if res.resname in ['PHE', 'DPHE', 'TYR', 'DTYR']:
            ring = res.atoms.select_atoms("name CG CD1 CD2 CE1 CE2 CZ")
            if len(ring) == 6:
                ring_atoms.append(ring)
    return ring_atoms

def calculate_normal_vector(ring):
    try:
        v1 = ring.positions[1] - ring.positions[0]
        v2 = ring.positions[2] - ring.positions[0]
        normal_vector = np.cross(v1, v2)
        norm = np.linalg.norm(normal_vector)
        return normal_vector / norm if norm > 0 else np.zeros(3)
    except Exception:
        return np.zeros(3)
