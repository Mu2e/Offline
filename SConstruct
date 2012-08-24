#!/usr/bin/env python
#
# Build a Mu2e base release or test release.
#
# $Id: SConstruct,v 1.39 2012/08/24 20:19:31 gandr Exp $
# $Author: gandr $
# $Date: 2012/08/24 20:19:31 $
#
# Original author Rob Kutschke.
#
import os
import sys

# Check that the release-specific setup has been run.
if not os.environ.has_key('MU2E_BASE_RELEASE'):
    sys.exit('You must setup a Mu2e base release before running scons.\nExiting.')

# Extract information from the shell environment.
art_inc       = os.environ['ART_INC']
art_lib       = os.environ['ART_LIB']
base          = os.environ['MU2E_BASE_RELEASE']
boost_lib     = os.environ['BOOST_LIB']
boost_inc     = os.environ['BOOST_INC']
clhep_inc     = os.environ['CLHEP_INC']
clhep_lib     = os.environ['CLHEP_LIB_DIR']
cppunit_dir   = os.environ['CPPUNIT_DIR']
gccxml_dir    = os.environ['GCCXML_DIR']
heppdt_lib    = os.environ['HEPPDT_LIB']
heppdt_inc    = os.environ['HEPPDT_INC']
libsigcpp_inc = os.environ['LIBSIGCPP_INC']
libsigcpp_lib = os.environ['LIBSIGCPP_LIB']
root_inc      = os.environ['ROOT_INC']
root_sys      = os.environ['ROOTSYS']
fhicl_inc     = os.environ['FHICLCPP_INC']
fhicl_lib     = os.environ['FHICLCPP_LIB']
cpp0x_inc     = os.environ['CPP0X_INC']
cpp0x_lib     = os.environ['CPP0X_LIB']
mesfac_inc     = os.environ['MESSAGEFACILITY_INC']
mesfac_lib     = os.environ['MESSAGEFACILITY_LIB']
cetlib_inc     = os.environ['CETLIB_INC']
cetlib_lib     = os.environ['CETLIB_LIB']
xercesc_inc    = os.environ['XERCES_C_INC']
xercesc_root   = os.environ['XERCESCROOT']

# If we are working in a test release, extract more information from the environment.
if os.environ.has_key('MU2E_TEST_RELEASE'):
    testrelease          = os.environ['MU2E_TEST_RELEASE']
    cpppath_frag         = [ testrelease, testrelease + '/BaBar/include' ]
    libpath_frag         = [ testrelease+'/lib/' ]
else:
    cpppath_frag         = [ ]
    libpath_frag         = [ ]

# The link libraries needed when building the BaBar code.
babarlibs = [ 'BaBar_KalmanTrack', 'BaBar_DetectorModel',  'BaBar_TrkBase',    'BaBar_BField',
              'BaBar_TrajGeom',    'BaBar_BbrGeom',        'BaBar_difAlgebra', 'BaBar_ProbTools',
              'BaBar_BaBar',       'BaBar_CLHEP',          'BaBar_MatEnv' ]

# Define scons-local environment - it will be exported later.

osenv = {}
for var in [ 'LD_LIBRARY_PATH',  'GCC_FQ_DIR',  'PATH', 'PYTHONPATH',  'ROOTSYS' ]:
    if var in os.environ.keys():
        osenv[var] = os.environ[var]
        pass
    pass

env = Environment( CPPPATH=[ cpppath_frag,
                             base,
                             base+'/BaBar/include',
                             art_inc,
                             mesfac_inc,
                             fhicl_inc,
                             cetlib_inc,
                             cpp0x_inc,
                             boost_inc,
                             clhep_inc,
                             cppunit_dir+'/include',
                             heppdt_inc,
                             libsigcpp_inc+'/sigc++-2.0',
                             libsigcpp_lib+'/sigc++-2.0/include',
                             root_inc,
                             xercesc_inc,
                           ],
                   LIBPATH=[ libpath_frag,
                             base+'/lib',
                             art_lib,
                             mesfac_lib,
                             fhicl_lib,
                             cetlib_lib,
                             cpp0x_lib,
                             boost_lib,
                             clhep_lib,
                             cppunit_dir+'/lib',
                             heppdt_lib,
                             libsigcpp_lib,
                             root_sys+'/lib',
                             '/lib', '/usr/X11R6/lib',
                             xercesc_root+'/lib',
                           ],
                   ENV=osenv,
                   FORTRAN = 'gfortran',
                   BABARLIBS = [ babarlibs ]
                 )

# Define the rule for building dictionaries.
genreflex_flags = '--deep --fail_on_warnings --iocomments --capabilities=classes_ids.cc '\
                + '-D_REENTRANT -DGNU_SOURCE -DGNU_GCC '\
                + '-DPROJECT_NAME="mu2e" -DPROJECT_VERSION="development"'
aa="if   t1=`expr ${TARGET} : '\(.*\)_dict.cpp'`;then t2=$${t1}_map.cpp; t1=$${t1}_dict.cpp;"\
  +"elif t1=`expr ${TARGET} : '\(.*\)_map.cpp'`;then t2=$${t1}_map.cpp; t1=$${t1}_dict.cpp; fi;"\
  +"if genreflex $SOURCE -s ${SOURCE.srcdir}/classes_def.xml $_CPPINCFLAGS"\
  +" -o $$t1 "\
  +genreflex_flags\
  +"; then mv ${TARGET.dir}/classes_ids.cc $$t2; else rm -f $$t1; false; fi"

genreflex = Builder(action=aa)
env.Append(BUILDERS = {'DictionarySource' : genreflex})

# Set compile and link flags.
SetOption('warn', 'no-fortran-cxx-mix')
env.MergeFlags('-g')
env.MergeFlags('-O3')
env.MergeFlags('-fno-omit-frame-pointer')
env.MergeFlags('-DNDEBUG')
env.MergeFlags('-rdynamic')
env.MergeFlags('-Wall')

# Extract gcc version.  Some libraries have this version embedded in their names.
ff = os.popen('g++ --version'); ll = ff.readline(); ff.close()
gcc_version = ll[10:13]
env.gcc_ver=gcc_version.replace('.','')

# Make the modified environment visible to all of the SConscript files
Export('env')

# Walk the directory tree to locate all SConscript files.
ss=[]
for root,dirs,files in os.walk('.'):
    for file in files:
        if file == 'SConscript': ss.append('%s/%s'%(root[2:],file))
        pass
    pass

# If the splines package is absent, skip building of the figure of merit tool.
if not os.environ.has_key('SPLINES_DIR'):
    if os.path.exists('FigureOfMerit/src/SConscript'):
        ss.remove('FigureOfMerit/src/SConscript')

# If the Dch part of BaBar package is not present, skip building KalmanTestsI.
if not(os.path.exists('BaBar/Dch')):	 
    ss.remove('KalmanTestsI/src/SConscript')
    print 'Dch part of the BaBar package is absent. Will not build KalmanTestsI.'

# Tell scons to operate on all of the SConscript files found in the previous steps.
env.SConscript(ss)

# This tells emacs to view this file in python mode.
# Local Variables:
# mode:python
# End:
