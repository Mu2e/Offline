# Ryunosuke O'Neil, 2020
# @ryuwd on GitHub
# millepede ---> proditions utility
# Parse millepede.res files and produce proditions text files


import io, sys

class AlignmentTable:
    def __init__(self, table, class_id, objects, ndof):
        self.table = table
        self.nobjects = objects
        self.ndof = ndof
        self.classid = class_id

        self.constants = []
        self.errors = []

        for _ in range(self.nobjects):
            self.constants.append([0.0]*ndof)
            self.errors.append([0.0]*ndof)

    def table_name(self):
        return self.table 
    
    def n_objects(self):
        return self.nobjects
    
    def dof_per_obj(self):
        return self.ndof
    
    def mplabel(self, id, dof):
        return 10000 * self.classid + 10 * id + dof
    
    def setConstant(self, id, dof, value):
        self.constants[id][dof] = float(value)
    
    def setError(self, id, dof, error):
        self.errors[id][dof] = float(error)

    def to_proditions_table(self):
        lines = []
        for i, constants in enumerate(self.constants):
            lines.append('%d,' % i + ','.join([str(i) for i in constants]))
        return """TABLE {name}
{csv}

""".format(name=self.table,
            csv='\n'.join(lines))

class AlignmentConstants:
    # interface to alignment constants conditions
    # and millepede results

    section_keyword = "TABLE"
    
    def __init__(self):
        self.tables = {'TrkAlignTracker': AlignmentTable('TrkAlignTracker', 0, 1, 6), 
                'TrkAlignPlane': AlignmentTable('TrkAlignPlane', 1, 36, 6), 
                'TrkAlignPanel': AlignmentTable('TrkAlignPanel', 2, 216, 6)}

    def read_db_file(self, input_file):
        section_name = ''
        with open(input_file, 'r') as f:
            for line in f.readlines():
                line = line.strip()

                if line.startswith('#'):
                    continue

                if len(line) == 0:
                    continue

                if '#' in line:
                    line = line.split('#')[0]

                if line.startswith(self.section_keyword):
                    section_name = line.split()[-1]
                    continue

                if section_name in self.tables:
                    splits = line.split(',')
                    if len(splits[1:]) == self.tables[section_name].dof_per_obj():
                        for i, val in enumerate(splits[1:]):
                            self.tables[section_name].setConstant(splits[0], i, val)


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
                cols = line.split()
                
                label, p, _ = cols[:3]
                p = float(p)
                perr = 0.0

                if len(cols) > 3:
                    perr = cols[-1]

                obj_type, id, dof =  self.parse_label(label)

                for table in self.tables.keys():
                    if obj_type == self.tables[table].classid:
                        self.tables[table].constants[id][dof]= p
                        self.tables[table].errors[id][dof] = perr
                        break

    def export_table(self):
        with io.StringIO() as f:
            for table in self.tables.values():
                f.write(table.to_proditions_table())
            return f.getvalue()


if __name__ == '__main__':
    input_file = 'millepede.res'

    if len(sys.argv) > 1:
        input_file = sys.argv[1]

    consts = AlignmentConstants()
    consts.read_mp_file(input_file)
    print (consts.export_table())
