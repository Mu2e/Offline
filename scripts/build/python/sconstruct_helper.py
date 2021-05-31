#
# Define some helper functions which provide text such as build options
# and library lists to be used in SConstruct.  Also there are a few functions
# that perform little tasks - put here to keep SConstruct more readable.
#

from glob import glob
import os, re, string

import sys
import subprocess

# Check that some of the required environment variables have been set
# and derive and check other pieces of the environment
# return a dictionary with mu2eOpts
def mu2eEnvironment():
    mu2eOpts = {}
    if 'MU2E_BASE_RELEASE' not in os.environ:
        raise Exception('You have not specified MU2E_BASE_RELEASE for this build')
    primaryBase = os.environ['MU2E_BASE_RELEASE']
    mu2eOpts["primaryBase"] = primaryBase

    if "MU2E_SATELLITE_RELEASE" in os.environ:
        mu2eOpts["satellite"] = True
        mu2eOpts["satelliteBase"] = os.environ["MU2E_SATELLITE_RELEASE"]
        base = mu2eOpts["satelliteBase"]
    else:
        mu2eOpts["satellite"] = False
        base = mu2eOpts["primaryBase"]

    # base is the local Offline area where output is written
    mu2eOpts["base"] = base
    mu2eOpts['libdir'] = base+'/lib'
    mu2eOpts['bindir'] = base+'/bin'
    mu2eOpts['tmpdir'] = base+'/tmp'
    mu2eOpts['gendir'] = base+'/gen'

    envopts = os.environ['MU2E_SETUP_BUILDOPTS'].strip()
    fsopts  = subprocess.check_output(primaryBase+"/buildopts",shell=True).strip().decode() # decode to convert byte string to text
    if envopts != fsopts:
        raise Exception("ERROR: Inconsistent build options: (MU2E_SETUP_BUILDOPTS vs ./buildopts)\n"
             +"Please source setup.sh after setting new options with buildopts.\n")

    # copy the buildopts to the dictionary
    mu2eOpts["buildopts"] = fsopts
    for line in fsopts.split():
        pp = line.split("=")
        mu2eOpts[pp[0]] = pp[1]  # e.g.  mu2eOpts["build"] = "prof"


    return mu2eOpts

# the list of root libraries
# This comes from: root-config --cflags --glibs
def rootLibs():
    return [ 'GenVector', 'Core', 'RIO', 'Net', 'Hist', 'MLP', 'Graf', 'Graf3d', 'Gpad', 'Tree',
             'Rint', 'Postscript', 'Matrix', 'Physics', 'MathCore', 'Thread', 'Gui', 'm', 'dl' ]


# the include path
def cppPath(mu2eOpts):
    path = [
        mu2eOpts["primaryBase"],
        os.environ['ART_INC'],
        os.environ['ART_ROOT_IO_INC'],
        os.environ['CANVAS_INC'],
        os.environ['BTRK_INC'],
        os.environ['KINKAL_INC'],
        os.environ['MESSAGEFACILITY_INC'],
        os.environ['FHICLCPP_INC'],
        os.environ['HEP_CONCURRENCY_INC'],
        os.environ['SQLITE_INC'],
        os.environ['CETLIB_INC'],
        os.environ['CETLIB_EXCEPT_INC'],
        os.environ['BOOST_INC'],
        os.environ['CLHEP_INC'],
        os.environ['HEPPDT_INC'],
        os.environ['ROOT_INC'],
        os.environ['XERCES_C_INC'],
        os.environ['TBB_INC'],
        os.environ['MU2E_ARTDAQ_CORE_INC'],
        os.environ['ARTDAQ_CORE_INC'],
        os.environ['PCIE_LINUX_KERNEL_MODULE_INC'],
        os.environ['TRACE_INC'],
        os.environ['GSL_INC'],
        os.environ['POSTGRESQL_INC']
        ]

    if mu2eOpts['satellite']:
        path = [ mu2eOpts['satelliteBase'] ] + path

    return path

# the ld_link_library path
def libPath(mu2eOpts):
    path = [
        mu2eOpts['primaryBase']+'/lib',
        os.environ['ART_LIB'],
        os.environ['ART_ROOT_IO_LIB'],
        os.environ['CANVAS_LIB'],
        os.environ['BTRK_LIB'],
        os.environ['KINKAL_LIB'],
        os.environ['MU2E_ARTDAQ_CORE_LIB'],
        os.environ['ARTDAQ_CORE_LIB'],
        os.environ['PCIE_LINUX_KERNEL_MODULE_LIB'],
        os.environ['MESSAGEFACILITY_LIB'],
        os.environ['HEP_CONCURRENCY_LIB'],
        os.environ['FHICLCPP_LIB'],
        os.environ['SQLITE_LIB'],
        os.environ['CETLIB_LIB'],
        os.environ['CETLIB_EXCEPT_LIB'],
        os.environ['BOOST_LIB'],
        os.environ['CLHEP_LIB_DIR'],
        os.environ['HEPPDT_LIB'],
        os.environ['ROOTSYS']+'/lib',
        os.environ['XERCESCROOT']+'/lib',
        os.environ['TBB_LIB'],
        os.environ['GSL_LIB'],
        os.environ['POSTGRESQL_LIBRARIES']
        ]

    if mu2eOpts['satellite']:
        path = [ mu2eOpts['satelliteBase']+'/lib' ] + path

    return path

# Define the compiler and linker options.
# These are given to scons using its Evironment.MergeFlags call.
def mergeFlags(mu2eOpts):
    build = mu2eOpts['build']
    flags = ['-std=c++17','-Wall','-Wno-unused-local-typedefs','-g',
             '-Werror','-Wl,--no-undefined','-gdwarf-2', '-Wl,--as-needed',
             '-Werror=return-type','-Winit-self','-Woverloaded-virtual', '-DBOOST_BIND_GLOBAL_PLACEHOLDERS' ]
    if build == 'prof':
        flags = flags + [ '-O3', '-fno-omit-frame-pointer', '-DNDEBUG' ]
    elif build == 'debug':
        flags = flags + [ '-O0' ]
    return flags


# Prepare some shell environmentals in a form to be pushed
# into the scons environment.
def exportedOSEnvironment():
    osenv = {}
    for var in [ 'LD_LIBRARY_PATH',  'GCC_FQ_DIR',  'PATH', 'PYTHONPATH',
                 'ROOTSYS', 'PYTHON_ROOT', 'PYTHON_DIR', 'SQLITE_FQ_DIR' ]:
        if var in os.environ.keys():
            osenv[var] = os.environ[var]
    return osenv

# list of BaBar libs
def BaBarLibs():
    return [ 'BTrk_KalmanTrack', 'BTrk_DetectorModel', 'BTrk_TrkBase',
             'BTrk_BField','BTrk_BbrGeom', 'BTrk_difAlgebra',
             'BTrk_ProbTools','BTrk_BaBar', 'BTrk_MatEnv' ]

# Walk the directory tree to locate all SConscript files.
def sconscriptList(mu2eOpts):
    ss=[]
    ss_append = ss.append
    for root, _, files in os.walk('.'):
        if 'SConscript' in files:
            ss_append(os.path.join(root[2:], 'SConscript'))

    return ss

# Make sure the build directories are created
def makeSubDirs(mu2eOpts):
    for mdir in [mu2eOpts[d] for d in ['libdir','bindir','tmpdir', 'gendir']]:
        os.makedirs(mdir, exist_ok=True)


# with -c, scons will remove all dependant files it knows about
# but when a source file is deleted:
# - the .os file is harmless since it will be ignored
# - the dict and lib now contain extra objects but scons can't
#     know that, so explicitly delete them here
def extraCleanup():
    for top, dirs, files in os.walk("./lib"):
        for name in files:
            ff =  os.path.join(top, name)
            print("removing file ", ff)
            os.unlink (ff)

    for top, dirs, files in os.walk("./tmp"):
        for name in files:
            ff =  os.path.join(top, name)
            print("removing file ", ff)
            os.unlink (ff)

    for top, dirs, files in os.walk("./gen"):
        for name in files:
            ff =  os.path.join(top, name)
            print("removing file ", ff)
            os.unlink (ff)
