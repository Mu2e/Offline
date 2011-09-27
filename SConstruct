# Build a Mu2e base release or test release.
#
# $Id: SConstruct,v 1.25 2011/09/27 21:48:31 mu2ecvs Exp $
# $Author: mu2ecvs $
# $Date: 2011/09/27 21:48:31 $
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
python_dir    = os.environ['PYTHON_DIR']
root_inc      = os.environ['ROOT_INC']
root_sys      = os.environ['ROOTSYS']
scons_dir     = os.environ['SCONS_DIR']
fhicl_inc     = os.environ['FHICLCPP_INC']
fhicl_lib     = os.environ['FHICLCPP_LIB']
cpp0x_inc     = os.environ['CPP0X_INC']
cpp0x_lib     = os.environ['CPP0X_LIB']
mesfac_inc     = os.environ['MESSAGEFACILITY_INC']
mesfac_lib     = os.environ['MESSAGEFACILITY_LIB']
cetlib_inc     = os.environ['CETLIB_INC']
cetlib_lib     = os.environ['CETLIB_LIB']

# If we are working in a test release, extract more information from the environment.
if os.environ.has_key('MU2E_TEST_RELEASE'):
    testrelease          = os.environ['MU2E_TEST_RELEASE']
    cpppath_frag         = [ testrelease, testrelease + '/BaBar/include' ]
    libpath_frag         = [ testrelease+'/lib/' ]
else:
    cpppath_frag         = [ ]
    libpath_frag         = [ ]

# Define scons-local environment - it will be exported later.
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
                             python_dir+'/include',
                             root_inc,
                             scons_dir+'/include',
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
                             python_dir+'/lib',
                             root_sys+'/lib',
                             scons_dir+'/lib',
                             '/lib', '/usr/X11R6/lib',
                           ],
                   ENV={
                         'LD_LIBRARY_PATH': os.environ['LD_LIBRARY_PATH'],
                         'GCC_FQ_DIR' : os.environ['GCC_FQ_DIR'], # For GCCXML
                         'PATH' : os.environ['PATH'], 
                         'PYTHONPATH' : os.environ['PYTHONPATH'],
                         'ROOTSYS' : os.environ['ROOTSYS']
                       },
                   FORTRAN = 'gfortran'
                 )

# Define the rule for building dictionaries.
genreflex_flags = '--deep --fail_on_warnings  --capabilities=classes_ids.cc '\
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
env.MergeFlags('-O0')
env.MergeFlags('-rdynamic')
env.MergeFlags('-Wall')

# Extract gcc version.  Some libraries have this version embedded in their names.
ff = os.popen('g++ --version'); ll = ff.readline(); ff.close()
gcc_version = ll[10:13]
env.gcc_ver=gcc_version.replace('.','')

# Make the modified environment visible to all of the SConscript files
Export('env')

# Find all Package/src/SConscript files
ss=[]
for root,dirs,files in os.walk('.'):
    for file in files:
        if file == 'SConscript': ss.append('%s/%s'%(root[2:],file))
        pass
    pass

# If the BaBar package is not present, do not make packages that depend on it.
# This needs to be maintained by hand.
if not(os.path.exists('BaBar/BaBar/src/SConscript')):
    ss.remove('KalmanTests/src/SConscript')
    ss.remove('TrkPatRec/src/SConscript')
    print 'BaBar package is absent. Will not build packages that depend on it.'
else:
#  Remove Dch code for now: not used by any officlal package
    for root,dirs,files in os.walk('BaBar/Dch'):
	for file in files:
	    if file == 'SConscript': ss.remove('%s/%s'%(root[0:],file))
	    #print '%s/%s'%(root[0:],file)
	    pass
	pass

#
# Tell scons to operate on all of the SConscript files found in the previous step.
env.SConscript(ss)

# This tells emacs to view this file in python mode.
# Local Variables:
# mode:python
# End:
