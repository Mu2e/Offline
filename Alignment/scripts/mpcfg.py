
import sys
import argparse
import io

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

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


def get_params(o_cls, param_inputs={}, fix=False, rotations=True):
    lines = ["Parameter"]
    dofs = DOF_TRANSLATIONS
    if rotations:
        dofs += DOF_ROTATIONS

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
        convergence_criteria='0.01'
    )

    write_output(args, buf)

# TODO: generation of constraints files with arg-sensitive options
# TODO: generation of parameter fixing files
def generate_constraint_cfg(args):
    constraints = get_constraints(1, rotations=(not args.no_rotations))

    if not args.no_panels:
        constraints += get_constraints(2, rotations=(not args.no_rotations))

    write_output(args, constraints)

def generate_param_cfg(args):
    parameters = ''

    # TODO: allow fine tuning of params (fixing, etc) using opts
    parameters += get_params(1, rotations=(not args.no_rotations))

    if not args.no_panels:
        parameters += get_params(2, rotations=(not args.no_rotations))

    write_output(args, parameters)

class AlignmentConstants:
    # interface to alignment constants conditions

    section_keyword = "TABLE"
    sections = []


    plane_constants = {}

    panel_constants = {}

    def __init__(self, input_filename):
        self._input_file = input_filename

        self.section_handlers = {
            'TrkAlignTracker':(lambda x: None, self.write_tracker_constants),
            'TrkAlignPlane': (self.parse_plane_constants, self.write_plane_constants),
            'TrkAlignPanel': (self.parse_panel_constants, self.write_panel_constants),
        }

        self._read_db_file()


    def _read_db_file(self):
        section_lines = []
        section_name = ''
        with open(self._input_file, 'r') as f:
            for line in f.readlines():
                line = line.strip()

                if line.startswith('#'):
                    continue

                if len(line) == 0:
                    continue

                if line.startswith(self.section_keyword):
                    if len(section_lines) > 0:
                        self.sections.append(section_lines)

                    section_lines = []
                    section_name = line.split()[-1]
                    continue

                if section_name in self.section_handlers:
                    self.section_handlers[section_name][0](line)

    def parse_plane_constants(self, line):
        row = [i.strip() for i in line.split(',')]

        id = int(row[0])

        if id >= 0 and id < 36:
            dx,dy,dz,a,b,g = [float(i) for i in row[1:]]
            self.plane_constants[id] = [dx,dy,dz,a,b,g]

    def write_plane_constants(self, fi):
        fi.write('%s %s\n' % (self.section_keyword, 'TrkAlignPlane'))
        for id, row in self.plane_constants.items():
            dofrow = ','.join(['%.8f' % (dof) for dof in row])
            fi.write('%d,%s\n' % (id,dofrow))
        fi.write('\n')

    def parse_panel_constants(self, line):

        # FIXME! not DRY-adhering
        row = [i.strip() for i in line.split(',')]

        id = int(row[0])

        if id >= 0 and id < 216:
            dx,dy,dz,a,b,g = [float(i) for i in row[1:]]
            self.panel_constants[id] = [dx,dy,dz,a,b,g]

    def write_panel_constants(self, fi):
        fi.write('%s %s\n' % (self.section_keyword, 'TrkAlignPanel'))
        for id, row in self.panel_constants.items():
            dofrow = ','.join(['%.8f' % (dof) for dof in row])
            fi.write('%d,%s\n' % (id,dofrow))
        fi.write('\n')


    def write_tracker_constants(self,fi):
        fi.write("""
TABLE TrkAlignTracker
0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
""")


    def parse_label(self, strlabel):
        # returns OBJECT TYPE (plane/panel), OBJECT ID and DOF index
        obj_type = strlabel[0]
        obj_id = strlabel[1:4]
        obj_dof = strlabel[4]

        return int(obj_type), int(obj_id), int(obj_dof)


    def apply_mp_result(self, inputfile):
        with open(inputfile, 'r') as f:
            for line in f.readlines():
                line = line.strip()
                if 'Parameter' in line:
                    continue
                label, p, _, dp, error = line.split()

                p = float(p)
                dp = float(dp)
                error = float(error)

                obj_type, id, dof =  self.parse_label(label)

                if obj_type == 1:
                    self.plane_constants[id][dof] -= dp
                elif obj_type == 2:
                    self.panel_constants[id][dof] -= dp

    def export_table(self):
        with io.StringIO() as f:
            for (_, writer) in self.section_handlers.values():
                writer(f)

            return f.getvalue()


    def export_df_planes(self):
        df = pd.DataFrame(columns=["label", "value"])

        for id, v in self.plane_constants.items():
            for i, dof in enumerate(v):
                df = df.append({
                    "label": get_label(1, id, i),
                    "value":  dof,
                    }, ignore_index=True)
        return df

def generate_alignment_constants(args):
    result_file = args.mp_res_file # whitespace separated columns
    last_align_constants = args.db_constants_file # csv-like format

    obj = AlignmentConstants(last_align_constants)
    obj.apply_mp_result(result_file)

    write_output(args, obj.export_table())


def alignment_const_compare(args):
    files = parse_input(args)

    rs = []
    lens=[]


    for f in files:
        obj = AlignmentConstants(f)

        df = obj.export_df_planes()
        lens.append(len(df))

        if len(rs) == 0:
            rs.append(np.arange(len(df)))

        r = [x + 0.25 for x in rs[-1]]
        rs.append(r)

        plt.bar(r, df['value'], width=0.25, label=f)
    plt.xlabel('group', fontweight='bold')
    plt.xticks([r + 0.25 for r in range(lens[0])], ['P%d' % i for i in range(36)])
    plt.legend()

    plt.show()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("operating_mode", type=str,
                        help="operating mode. can be: steer, constr, param")

    parser.add_argument("-o", "--output", type=str, default='-',
                        help="output file")

    parser.add_argument("-v", "--verbose", action="store_true")

    parser.add_argument("-i", "--input", type=str, default='-',
                        help="input file(s). Can be binary or txt configs. Default: '-' for stdin (NEWLINE separated)")

    parser.add_argument("--no-panels", action="store_true")
    parser.add_argument("--no-rotations", action="store_true")

    parser.add_argument("--mp-res-file", type=str,
                        default='millepede.res',
                        help="millepede.res file from PEDE (mergeres only!)")

    parser.add_argument("--db-constants-file", type=str,
                        help="alignment constants DB file for the last track run iteration (mergeres only!)")

    args = parser.parse_args()

    if args.operating_mode == 'steer':
        generate_steering(args)

    elif args.operating_mode == 'constr':
        generate_constraint_cfg(args)

    elif args.operating_mode == 'param':
        generate_param_cfg(args)

    elif args.operating_mode == 'mergeres':
        if args.db_constants_file is None:
            parser.print_help()
            return

        if args.mp_res_file is None:
            parser.print_help()
            return

        generate_alignment_constants(args)

    elif args.operating_mode == 'constcomp':
        alignment_const_compare(args)

    else:
        parser.print_help()

if __name__ == '__main__':
    main()