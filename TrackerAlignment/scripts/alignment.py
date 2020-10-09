from sympy import Symbol, Matrix, diff, sqrt, atan2, cos, sin

def unit_vector(v):
    tot2 = v.dot(v)
    return v/sqrt(tot2)

# driftvel comes from driftInstantSpeed
def D2T(dca, driftvel):
    return dca / driftvel


def HitAmbiguity(intercept, tdir, straw_mp, straw_dir):
    sep = intercept - straw_mp
    perp = unit_vector(tdir.cross(straw_dir))
    dperp = perp.dot(sep)

    return dperp
    # to get the right sign: dperp > 0 ? -1 : 1

class TwoLinePCA:
    def __init__(self, p1, t1, p2, t2):
        t1 = unit_vector(t1)
        t2 = unit_vector(t2)

        c = t1.dot(t2)

        sinsq = 1.0 - c*c
        _delta = p1 - p2
        ddotT1 = _delta.dot(t1)
        ddotT2 = _delta.dot(t2)

        _s1 = (ddotT2*c-ddotT1)/sinsq
        _s2 = -(ddotT1*c-ddotT2)/sinsq

        _pca1 = p1 + t1 * _s1
        _pca2 = p2 + t2 * _s2

        pdiff = _pca1 - _pca2
        
        self._s1 = _s1
        self._s2 = _s2
        self._p1 = _pca1
        self._p2 = _pca2
        self._dca = sqrt(pdiff.dot(pdiff))

    def dca(self):
        return self._dca
    
    def s1(self):
        return self._s1
    
    def s2(self):
        return self._s2
    
    def point1(self):
        return self._p1
    
    def point2(self):
        return self._p2

def colvec_perp(matrix):
    return sqrt(matrix[0]*matrix[0] + matrix[1]*matrix[1])


def exact_rotation_matrix(a, b, g):
    R_x = Matrix([[1, 0, 0],
                  [0, cos(a), -sin(a)],
                  [0, sin(a), cos(a)]])

    R_y = Matrix([[cos(b), 0, sin(b)],
                  [0,      1, 0],
                  [-sin(b),0, cos(b)]])

    R_z = Matrix([[cos(g), -sin(g), 0],
                  [sin(g), cos(g),  0],
                  [0,      0,       1]])

    return R_z*(R_y*R_x)

# Mu2eUtilities.
class HepTransform:
    def __init__(self, trl, rot, mat=False):
        self._trl = trl

        if mat:
            self._rotmat = rot
            return

        a, b, g = rot

        self._rotmat = exact_rotation_matrix(a, b, g)

    def trl(self):
        return self._trl

    def rotmat(self):
        return self._rotmat

    def combine(self, b):
        v = self._trl + self._rotmat * b._trl
        r = self._rotmat * b._rotmat
        return HepTransform(v, r, mat=True)

    def transform(self, vector):
        return self._rotmat * vector + self._trl

# based on AlignedTrackerMaker in TrackerConditions.
def nested_transform_alignment(wire_pos, wire_dir,
                               plane_origin, panel_straw0mp,
                               plane_trl, plane_rot, panel_trl, panel_rot):

    align_plane = HepTransform(plane_trl, plane_rot)
    align_panel = HepTransform(panel_trl, panel_rot)

    plane_to_tracker = HepTransform(
        Matrix([0, 0.0, plane_origin[2]]), [0.0, 0.0, 0.0])

    dv = panel_straw0mp - Matrix([0, 0.0, plane_origin[2]])

    rz = 0
    if not (dv[0] == 0.0 and dv[1] == 0.0):
        rz = atan2(dv[1], dv[0])

    panel_to_plane = HepTransform(dv, [0.0, 0.0, rz])

    transform = plane_to_tracker.combine(
        align_plane.combine(panel_to_plane.combine(align_panel)))

    dx = colvec_perp(wire_pos) - colvec_perp(panel_straw0mp)
    dz = (wire_pos - panel_straw0mp)[2]

    straw_to_panel = Matrix([dx, 0.0, dz])
    straw_dir = Matrix([0.0, 1.0, 0.0])

    wire_pos_a = transform.transform(straw_to_panel)
    wire_dir_a = transform.rotmat() * straw_dir

    return wire_pos_a, wire_dir_a

