# Ryunosuke O'Neil, 2019
# Symbolic derivation of DOCA partial derivatives with respect to alignment and track params

# TODO: support Kalman tracks (KinKal?)
# TODO: cleanup
# TODO: unit testing
# TODO: DOCAtoTOCA should be a replica of D2T which properly accounts for nonlinear drift effects

import sympy
from sympy import Symbol, Matrix, diff, sqrt, atan2, cos, sin
from sympy.physics.vector import ReferenceFrame
from sympy.vector import matrix_to_vector, CoordSys3D
from sympy.functions import sign

from sympy.vector.orienters import BodyOrienter
from sympy.simplify.cse_main import cse
from sympy.utilities.codegen import codegen
from sympy.printing import ccode
from sympy.utilities.iterables import numbered_symbols
from sympy.functions import Abs
from sympy.matrices.dense import matrix_multiply_elementwise


def unit_vector(v):
    tot2 = v.dot(v)
    return v/sqrt(tot2)


def DOCAToTOCA(dca):
    return dca / 0.0625


def DOCA(p1, t1, p2, t2):
    t1 = unit_vector(t1)
    t2 = unit_vector(t2)
    # t2 should already be a unit vector

    c = t1.dot(t2)

    sinsq = 1.0 - c*c
    _delta = p1 - p2
    ddotT1 = _delta.dot(t1)
    ddotT2 = _delta.dot(t2)

    _s1 = (ddotT2*c-ddotT1)/sinsq
    _s2 = -(ddotT1*c-ddotT2)/sinsq

    _pca1 = p1 + t1 * _s1
    _pca2 = p2 + t2 * _s2

    ___diff = _pca1 - _pca2

    dca = sqrt(___diff.dot(___diff))

    return sympy.Piecewise((dca, _s2 > 0), (-1.0*dca, True))


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


class HepTransform:
    def __init__(self, trl, rot, mat=False):
        self._trl = trl

        if mat:
            self._rotmat = rot
            return

        a, b, g = rot

        # This is an approximation.
        # self._rotmat = Matrix([[1, g, b],
        #                        [-g, 1, a],
        #                        [b, -a, 1]])

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


def generate_expressions(approximate=False, remove_globalparam_dependence=True, time_domain=True):
    # define symbols for alignment and track parameters

    # plane alignment

    # translation vector
    dx = Symbol('plane_dx', real=True)
    dy = Symbol('plane_dy', real=True)
    dz = Symbol('plane_dz', real=True)
    plane_trl = Matrix([dx, dy, dz])

    # rotation angles
    a = Symbol('plane_a', real=True)
    b = Symbol('plane_b', real=True)
    g = Symbol('plane_g', real=True)
    plane_rot = (a, b, g)

    # panel alignment
    panel_a = Symbol('panel_a', real=True)
    panel_b = Symbol('panel_b', real=True)
    panel_g = Symbol('panel_g', real=True)
    panel_rot = (panel_a, panel_b, panel_g)

    # translation vector
    panel_dx = Symbol('panel_dx', real=True)
    panel_dy = Symbol('panel_dy', real=True)
    panel_dz = Symbol('panel_dz', real=True)
    panel_trl = Matrix([panel_dx, panel_dy, panel_dz])

    # track parametrisation
    a0 = Symbol('a0', real=True)
    b0 = Symbol('b0', real=True)
    track_pos = Matrix([a0, 0, b0])

    a1 = Symbol('a1', real=True)
    b1 = Symbol('b1', real=True)
    track_dir = Matrix([a1, -1, b1])

    t0 = Symbol('t0')

    # wire position (midpoint) and direction
    wx = Symbol('wire_x', real=True)
    wy = Symbol('wire_y', real=True)
    wz = Symbol('wire_z', real=True)
    wire_pos = Matrix([wx, wy, wz])

    wwx = Symbol('wdir_x', real=True)
    wwy = Symbol('wdir_y', real=True)
    wwz = Symbol('wdir_z', real=True)
    wire_dir = Matrix([wwx, wwy, wwz])

    # origin of the plane being aligned
    ppx = Symbol('plane_x', real=True)
    ppy = Symbol('plane_y', real=True)
    ppz = Symbol('plane_z', real=True)
    plane_origin = Matrix([ppx, ppy, ppz])

    # origin of the panel being aligned
    panel_x = Symbol('panel_straw0x', real=True)
    panel_y = Symbol('panel_straw0y', real=True)
    panel_z = Symbol('panel_straw0z', real=True)
    panel_straw0mp = Matrix([panel_x, panel_y, panel_z])

    local_params = [a0, b0, a1, b1, t0]
    global_params = [dx, dy, dz, a, b, g]
    global_params += [panel_dx, panel_dy, panel_dz, panel_a, panel_b, panel_g]

    wire_params = [wx, wy, wz, wwx, wwy, wwz]
    plane_position = [ppx, ppy, ppz]
    panel_position = [panel_x, panel_y, panel_z]

    all_params = local_params + global_params + \
        wire_params + plane_position + panel_position

    param_dict = {
        'all': all_params,
        'local': local_params,
        'global': global_params,
        'wire': wire_params,
        'plane_pos': plane_position,
        'panel_pos': panel_position
    }

    # choose method to align vectors with
    alignment_func = nested_transform_alignment

    # recalculate wire position and rotation according to alignment parameters
    aligned_wpos, aligned_wdir = alignment_func(
        wire_pos, wire_dir,
        plane_origin, panel_straw0mp,
        plane_trl, plane_rot,
        panel_trl, panel_rot)

    # this is the residual expression (excluding terms with no dependence on local
    # and global parameters )
    aligned_doca = aligned_doca_to_diff = DOCA(aligned_wpos, aligned_wdir, track_pos, track_dir)

    if time_domain:
        # we convert the DOCA to a TOCA and add T0 (since it is a local param)
        # the hit time has no explicit track parameter dependence
        aligned_doca_to_diff = DOCAToTOCA(aligned_doca) - t0

    # now generate optimised C code to calculate each deriv
    if remove_globalparam_dependence:
        param_dict['all'] = local_params + \
            wire_params + plane_position + panel_position

    expressions = []
    for parameter in local_params + global_params:
        # calculate derivative symbolically for each local and global parameter then generate code for the function
        pdev = diff(aligned_doca_to_diff, parameter)

        if remove_globalparam_dependence:
            # since these derivatives are evaluated with alignment dofs kept to zero
            # we substitute zero for
            # a, b, g, dx, dy, dz in our final expressions
            pdev = pdev.subs({
                dx: 0, dy: 0, dz: 0,
                a: 0, b: 0, g: 0,
                panel_dx: 0,
                panel_dy: 0,
                panel_dz: 0,
                panel_a: 0,
                panel_b: 0,
                panel_g: 0
            })
        expressions.append(pdev)

    # if remove_globalparam_dependence:
    #     aligned_doca = aligned_doca.subs({
    #             dx: 0, dy: 0, dz: 0,
    #             a: 0, b: 0, g: 0,
    #             panel_dx: 0,
    #             panel_dy: 0,
    #             panel_dz: 0,
    #             panel_a: 0,
    #             panel_b: 0,
    #             panel_g: 0
    #         })

    nominal_doca = DOCA(wire_pos, wire_dir, track_pos, track_dir)

    if VALIDATE:

        print(aligned_wpos)

        aligned_wpos = aligned_wpos.subs({
            dx: 1.50269,
            dy: -1.94275,
            dz: 3.27363,

            a: 0,
            b: 0,
            g: 0,

            panel_dx: 0,
            panel_dy: 0,
            panel_dz: 0,
            panel_a: 0,
            panel_b: 0,
            panel_g: 0,

            wx: 367.052,  # wire pos before (straw 0 plane 0 panel 0)
            wy: 98.3512,
            wz: -1490.37,

            panel_x: 368.561,  # panel straw0mp
            panel_y: 98.7556,
            panel_z: -1493.08,

            ppx: 0,  # plane origin
            ppy: 0,
            ppz: -1504.35
        })

        print(aligned_wpos)

        print('should be (368.555,96.4085,-1487.1)')

    return expressions, param_dict, aligned_doca, aligned_wpos, aligned_wdir


VALIDATE = False

c_template = """

# include "TrackerAlignment/inc/AlignmentDerivatives.hh"
# include <math.h>
# include <vector>

%s

"""

h_template = """
# ifndef RIGIDBODYDOCADERIV_H
# define RIGIDBODYDOCADERIV_H
# include <vector>
%s

# endif

"""

fn_template = "%s %s(%s);"


def cseexpr_to_ccode(symname, symfunc, symbolslist):
    tmpsyms = numbered_symbols("R")
    symbols, simple = cse(symfunc, symbols=tmpsyms)

    code = "double %s(%s)\n" % (str(symname), ", ".join(
        "double %s" % x for x in symbolslist))
    code += "{\n"
    for s in symbols:
        code += "    double %s = %s;\n" % (ccode(s[0]), ccode(s[1]))
    code += "    double result = %s;\n" % ccode(simple[0])
    code += "    return result;\n"
    code += "}\n"

    return code


def build_ccode_function(return_type, fn_name, symbolslist, fn_body):
    args = 'double ' + ', double '.join([symb.name for symb in symbolslist])
    code = """{return_type} {fn_name}({arg_list})
{{
        {body}
}}""".format(
        return_type=return_type,
        fn_name=fn_name,
        arg_list=args,
        body=fn_body
    )

    function_header_code = fn_template % (return_type, fn_name, args)

    return function_header_code, code


def generate_code_function(name, return_type, expr, symbols):
    args = 'double %s' % (', double '.join([p.name for p in symbols]))

    function_header_code = fn_template % (return_type, name, args)
    function_code = cseexpr_to_ccode(name, expr, symbols)

    return function_header_code, function_code


def main():
    function_prefix = "CosmicTrack_DCA"

    exprs, params, aligned_doca, aligned_wpos, aligned_wdir = generate_expressions()

    if VALIDATE:
        return

    lgparams = params['local'] + params['global']

    generated_code = []

    # generate code for DOCA calculation ( no global parameter dependence )
    generated_code.append(generate_code_function(
        function_prefix, 'double', aligned_doca, params['all'] + params['global']))

    # generate code for alignment function
    generated_code.append(generate_code_function(
        function_prefix+'alignpos_x', 'double', aligned_wpos[0], params['global'] + params['plane_pos'] + params['panel_pos'] + params['wire']))
    generated_code.append(generate_code_function(
        function_prefix+'alignpos_y', 'double', aligned_wpos[1], params['global'] + params['plane_pos'] + params['panel_pos'] + params['wire']))
    generated_code.append(generate_code_function(
        function_prefix+'alignpos_z', 'double', aligned_wpos[2], params['global'] + params['plane_pos'] + params['panel_pos'] + params['wire']))

    generated_code.append(generate_code_function(
        function_prefix+'aligndir_x', 'double', aligned_wdir[0], params['global'] + params['plane_pos'] + params['panel_pos'] + params['wire']))
    generated_code.append(generate_code_function(
        function_prefix+'aligndir_y', 'double', aligned_wdir[1], params['global'] + params['plane_pos'] + params['panel_pos'] + params['wire']))
    generated_code.append(generate_code_function(
        function_prefix+'aligndir_z', 'double', aligned_wdir[2], params['global'] + params['plane_pos'] + params['panel_pos'] + params['wire']))

    # generate code for all expressions
    generated_code += [generate_code_function('{}_Deriv_{}'.format(
        function_prefix, symb.name), 'double', expr, params['all']) for expr, symb in zip(exprs, lgparams)]

    # generate code to build arrays of calculated derivatives

    def generate_dcaderiv_arraybuilder_fn(params_type):
        code = "std::vector<double> result = {"
        code += ',\n'.join(['{fn}({args})'.format(
                fn=function_prefix + '_Deriv_' + symb.name,
                args=','.join([p.name for p in params['all']])
        ) for symb in params[params_type]
        ])
        code += '};\nreturn result;'
        return [build_ccode_function(
            'std::vector<double>', '{}_{}Deriv'.format(function_prefix, params_type.capitalize()), params['all'], code)]

    generated_code += generate_dcaderiv_arraybuilder_fn('local')
    generated_code += generate_dcaderiv_arraybuilder_fn('global')

    # pull everything together and write to file
    functions, code = zip(*generated_code)
    c_code = c_template % ('\n\n'.join(code))
    c_header = h_template % '\n\n'.join(functions)

    with open('src/AlignmentDerivatives.cc', 'w') as f:
        f.write(c_code)

    with open('inc/AlignmentDerivatives.hh', 'w') as f:
        f.write(c_header)


if __name__ == "__main__":
    main()
