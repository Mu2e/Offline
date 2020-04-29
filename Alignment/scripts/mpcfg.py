
import sys
import argparse

STEER_SKELETON = """Cfiles
{datafiles}

{constraints}

{parameters}

outlierdownweighting 4

method {method} {n_iterations} {convergence_criteria}

end
"""

def get_label(o_cls, planeno, dof):
    return (o_cls*10000+ 10*planeno+ dof)

DOF_TRANSLATIONS = [0,1,2]
DOF_ROTATIONS = [3,4,5]

OBJECTS ={
    1: range(0, 36),
    2: range(0, 216)
}



def get_constraints(o_cls, translations=True, rotations=True):
    lines = []
    dofs = []
    if translations: dofs += DOF_TRANSLATIONS
    if rotations: dofs += DOF_ROTATIONS

    for dof in DOF_TRANSLATIONS:
        lines.append("Constraint    0")
        for plane in OBJECTS[o_cls]:
            lines.append("%d    1" % get_label(o_cls, plane, dof))

    return '\n'.join(lines) + '\n'


def get_params(o_cls, param_inputs={}, fix=False):
    lines = ["Parameter"]
    dofs = DOF_TRANSLATIONS + DOF_ROTATIONS
    defps = 0
    if fix: defps = -0.5
    for dof in dofs:
        for plane in OBJECTS[o_cls]:
            lbl = get_label(o_cls, plane, dof)
            param_ins = (0, defps)
            if lbl in param_inputs: param_ins = param_inputs[lbl]
            lines.append("%d    %.2f    %.2f" % (lbl, *param_ins))

    return '\n'.join(lines) + '\n'


def parse_input(args):
    if args.input == '-':
        input_str = sys.stdin.read()
        return input_str.split()

    return args.input.split()

def write_output(args, buf):

    if args.output == '-':
        sys.stdout.write(buf)
    else:
        with open(args.output, 'w') as f:
            f.write(buf)

def generate_steering(args):
    files = parse_input(args)

    buf = STEER_SKELETON.format(
        datafiles='\n'.join(files),
        constraints='',
        parameters='',
        method='inversion',
        n_iterations='10',
        convergence_criteria='0.1'
    )

    write_output(args, buf)

# TODO: generation of constraints files with arg-sensitive options
# TODO: generation of parameter fixing files
def generate_constraint_cfg(args):
    constraints = get_constraints(1) + get_constraints(2)

    write_output(args, constraints)

def generate_param_cfg(args):
    parameters = ''

    opts = args.constr_dofs.split(',')

    # TODO: allow fine tuning of params (fixing, etc) using opts
    parameters += get_params(1)
    parameters += get_params(2, fix=True)

    write_output(args, parameters)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("operating_mode", type=str,
                        help="operating mode. can be: steer, constr, param")

    group = parser.add_mutually_exclusive_group()
    group.add_argument("--constr-dofs", type=str,
                        help='Constrain average rotation or translation to zero. comma separated list containing any or all of: all, pl-[trl/rot], pa-[trl/rot]')

    group.add_argument("--fix-dofs", type=str,
                        help='comma separated list containing any or all of: plane-[trl/rot], panel-[trl/rot]')

    parser.add_argument("-o", "--output", type=str, default='-',
                        help="output file")

    parser.add_argument("-v", "--verbose", action="store_true")

    parser.add_argument("-i", "--input", type=str, default='-',
                        help="input file(s). Can be binary or txt configs. Default: '-' for stdin (NEWLINE separated)")

    parser.add_argument("-g", "--gen-constr", action="store_true",
                        help="include constraints and parameter fixing information in the steering file")

    args = parser.parse_args()

    if args.operating_mode == 'steer':
        generate_steering(args)

    elif args.operating_mode == 'constr':
        generate_constraint_cfg(args)

    elif args.operating_mode == 'param':
        generate_param_cfg(args)

    else:
        parser.print_help()

if __name__ == '__main__':
    main()