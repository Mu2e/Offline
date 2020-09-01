
import sys, os
import pathlib
def scan_include(line):
    # scrape path
    mu2e_base = os.environ['MU2E_BASE_RELEASE']
    line = line.rsplit('//',1)[0]
    line = line.replace('#include', '').strip()
    line = line.replace('"', '').strip()
    line = line.replace('<', '').strip()
    line = line.replace('>', '').strip()

    abspath = os.path.join(mu2e_base, line)

    if os.path.exists(abspath) and os.path.isfile(abspath):
        return 'mu2e_%s' % pathlib.Path(line).parts[0]

    return None


def scan_file(fname):
    pkgs = []
    with open(fname, 'r') as f:
        for line in f:
            if line.startswith('#include'):
                pkg = scan_include(line.strip())

                if pkg is not None:
                    pkgs.append(pkg)
    return pkgs


def main():
    package = sys.argv[1]
    package_src = os.path.join(package, 'src/')

    files = [os.path.join(package_src, f) for f in os.listdir(package_src) if os.path.isfile(os.path.join(package_src, f))]

    packages = []
    for filename in files:
        if filename.endswith('.cc'):
            # scan includes
            packages += (scan_file(filename))

    print( 'libs = [')
    for package in set(packages):
        print ("    '%s'," % package)
    print( ']')

if __name__ == '__main__':
    main()