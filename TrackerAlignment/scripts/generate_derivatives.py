# Ryunosuke O'Neil, 2019
# Symbolic derivation of DOCA partial derivatives with respect to alignment and track params

# TODO: support Kalman tracks (KinKal?)
# TODO: cleanup
# TODO: unit testing

import os

import sympy
from sympy import Symbol, Matrix, diff, sqrt, atan2, cos, sin

from sympy.printing import ccode
from sympy.simplify.cse_main import cse
from sympy.utilities.codegen import codegen
from sympy.utilities.iterables import numbered_symbols


from alignment import nested_transform_alignment, unit_vector, TwoLinePCA, HitAmbiguity, D2T


def generate_expressions(approximate=False, remove_globalparam_dependence=False, time_domain=True):
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
    track_dir = unit_vector(Matrix([a1, -1, b1]))

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


    driftvel = Symbol('driftvel', real=True)

    local_params = [a0, b0, a1, b1, t0]
    global_params = [dx, dy, dz, a, b, g]
    global_params += [panel_dx, panel_dy, panel_dz, panel_a, panel_b, panel_g]

    wire_params = [wx, wy, wz, wwx, wwy, wwz]
    plane_position = [ppx, ppy, ppz]
    panel_position = [panel_x, panel_y, panel_z]

    all_params = local_params + global_params + \
        wire_params + plane_position + panel_position + [driftvel]

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

    # this is the predicted function expression 'f'
    pca = TwoLinePCA(track_pos, track_dir, aligned_wpos, aligned_wdir)
    ambig = HitAmbiguity(track_pos,track_dir,aligned_wpos,aligned_wdir)
    aligned_doca = diff_expr = sympy.Piecewise((-pca.dca(), ambig > 0), (pca.dca(), True)) 

    if time_domain:
        # FIXME! need to add time offset
        traj_time = (pca.point1() - track_pos).dot(track_dir)
        diff_expr = D2T(pca.dca(), driftvel) + t0 + traj_time #+ time_offset

    # now generate optimised C code to calculate each deriv
    if remove_globalparam_dependence:
        param_dict['all'] = local_params + \
            wire_params + plane_position + panel_position

    expressions = []
    for parameter in local_params + global_params:
        # calculate derivative symbolically for each local and global parameter then generate code for the function
        pdev = diff(diff_expr, parameter)
        expressions.append(pdev)

    return expressions, param_dict, aligned_doca, aligned_wpos, aligned_wdir


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




def fmt_function(return_type, fn_name, symbolslist, fn_body):
    args = 'double const& ' + ', double const& '.join([symb.name for symb in symbolslist])
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


def fn_from_expr(name, return_type, expr, symbols):
    def build_cse_fn(symname, symfunc, symbolslist):
        tmpsyms = numbered_symbols("R")
        symbols, simple = cse(symfunc, symbols=tmpsyms)

        code = "double %s(%s)\n" % (str(symname), ", ".join(
            "double const& %s" % x for x in symbolslist))
        code += "{\n"
        for s in symbols:
            code += "    double %s = %s;\n" % (ccode(s[0]), ccode(s[1]))
        code += "    double result = %s;\n" % ccode(simple[0])
        code += "    return result;\n"
        code += "}\n"

        return code

    args = 'double const& %s' % (', double const& '.join([p.name for p in symbols]))

    function_header_code = fn_template % (return_type, name, args)
    function_code = build_cse_fn(name, expr, symbols)

    return function_header_code, function_code


def main():
    function_prefix = "GaussianDriftFit_ResidDeriv_"

    if 'MU2E_BASE_RELEASE' not in os.environ:
        print ("Please source setup.sh for the Offline build you want to generate the code for.")
        return
    
    base_path = os.environ['MU2E_BASE_RELEASE']
    package_name = 'TrackerAlignment'

    print ("MU2E_BASE_RELEASE: %s" % base_path)
    print ("Package: %s" % package_name)

    base_path = os.path.join(base_path, package_name)

    print ("algebra...")
    exprs, params, aligned_doca, _, _ = generate_expressions()


    print ("code generation...")
    generated_code = []
    lgparams = params['local'] + params['global']

    # generate code for DOCA calculation ( no global parameter dependence )
    generated_code.append(fn_from_expr(
        function_prefix, 'double', aligned_doca, params['all']))

    # generate code for all expressions
    generated_code += [fn_from_expr('{}_Deriv_{}'.format(
        function_prefix, symb.name), 'double', expr, params['all']) for expr, symb in zip(exprs, lgparams)]

    # generate code to build arrays of calculated derivatives
    def generate_deriv_array_fn(params_type):
        code = "std::vector<double> result = {"
        code += ',\n'.join(['{fn}({args})'.format(
                fn=function_prefix + '_Deriv_' + symb.name,
                args=','.join([p.name for p in params['all']])
        ) for symb in params[params_type]
        ])
        code += '};\nreturn result;'
        return [fmt_function(
            'std::vector<double>', '{}_{}Deriv'.format(function_prefix, params_type.capitalize()), params['all'], code)]

    generated_code += generate_deriv_array_fn('local')
    generated_code += generate_deriv_array_fn('global')

    # pull everything together and write to file
    functions, code = zip(*generated_code)
    c_code = c_template % ('\n\n'.join(code))
    c_header = h_template % '\n\n'.join(functions)

    source_filename = os.path.join(base_path, 'src/AlignmentDerivatives.cc')
    header_filename = os.path.join(base_path, 'inc/AlignmentDerivatives.hh')

    with open(source_filename, 'w') as f:
        print ("writing %s.." % source_filename)
        f.write(c_code)

    with open(header_filename, 'w') as f:
        print ("writing %s.." % header_filename)
        f.write(c_header)

    print ("Wrote %d generated functions successfully. Please validate before using in production." % len(functions))

if __name__ == "__main__":
    main()
