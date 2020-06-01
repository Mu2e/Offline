import sys
# merge steering files into one where
# all input files (.bin) are included

def main():
    files=sys.argv[1:]
    file_lines=[]

    nlines = -1

    for filename in files:
        with open(filename, 'r') as f:
            lines = f.readlines()
            if nlines == -1:
                nlines = len(lines)
            if nlines != len(lines) or nlines < 0:
                exit(1)
            file_lines.append(lines)

    for lineno in range(nlines):
        lines = [file_lines[i][lineno] for i, _ in enumerate(files)]
        if len(lines) == 0:
            continue
        if '.bin' in lines[0]: # include all mille file inputs
            for line in lines:
                print (line.strip())
        else:
            pline = 0
            for line in lines:
                # assert all other config lines equal
                if pline != 0 and line != pline and '.txt' not in line:
                    exit(1)
                pline = line
            print(lines[0].strip())

if __name__ == '__main__':
    main()