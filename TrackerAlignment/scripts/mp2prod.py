import io

class AlignmentConstants:
    # interface to alignment constants conditions
    # and millepede results

    section_keyword = "TABLE"
    sections = []

    plane_constants = {}

    panel_constants = {}

    def __init__(self):
        self.section_handlers = {
            'TrkAlignTracker':(lambda x: None, self.write_tracker_constants),
            'TrkAlignPlane': (self.parse_plane_constants, self.write_plane_constants),
            'TrkAlignPanel': (self.parse_panel_constants, self.write_panel_constants),
        }


    def read_db_file(self, input_file):
        section_lines = []
        section_name = ''
        with open(input_file, 'r') as f:
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


    def read_mp_file(self, inputfile):
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
                    self.plane_constants[id][dof] = dp
                elif obj_type == 2:
                    self.panel_constants[id][dof] = dp

    def export_table(self):
        with io.StringIO() as f:
            for (_, writer) in self.section_handlers.values():
                writer(f)

            return f.getvalue()


if __name__ == '__main__':
    input_file = 'millepede.res'
    consts = AlignmentConstants()
    consts.read_mp_file(input_file)
    consts.export_table()
