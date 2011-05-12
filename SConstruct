# Build a Mu2e base release or test release.
#
# $Id: SConstruct,v 1.12 2011/05/12 15:50:46 kutschke Exp $
# $Author: kutschke $
# $Date: 2011/05/12 15:50:46 $
#
# Original author Rob Kutschke.
#
# Notes:
# 1) When one is not working in a base release, the variables testrelease*
#    are set to None.  So CPPPATH and LIBPATH will contain some entries
#    that evaluate to None.  There is some magic inside python that
#    skips the generation of -I and -L options for these files.  This
#    would also work had we set the variables to '' rather than to None.
#
import os
import sys

# Check that the site-local setup has been run.
if not os.environ.has_key('FRAMEWORK_DIR'):
    sys.exit('You must setup the externals before running scons.\nExiting.')

# Check that the release-specific setup has been run.
if not os.environ.has_key('MU2E_BASE_RELEASE'):
    sys.exit('You must setup a Mu2e base release before running scons.\nExiting.')

# Extract information from the shell environment.
framework     = os.environ['FRAMEWORK_DIR']
base          = os.environ['MU2E_BASE_RELEASE']
boost_lib     = os.environ['BOOST_LIB']
boost_inc     = os.environ['BOOST_INC']
clhep_inc     = os.environ['CLHEP_INC']
clhep_base    = os.environ['CLHEP_BASE']
cppunit_dir   = os.environ['CPPUNIT_DIR']
gccxml_dir    = os.environ['GCCXML_DIR']
heppdt_lib    = os.environ['HEPPDT_LIB']
heppdt_inc    = os.environ['HEPPDT_INC']
libsigcpp_inc = os.environ['LIBSIGCPP_INC']
libsigcpp_lib = os.environ['LIBSIGCPP_LIB']
python_dir    = os.environ['PYTHON_DIR']
root_dir      = os.environ['ROOT_DIR']
scons_dir     = os.environ['SCONS_DIR']

# If we are working in a test release, extract more information from the environment.
# See note 1.
if os.environ.has_key('MU2E_TEST_RELEASE'):
    testrelease          = os.environ['MU2E_TEST_RELEASE']
    testrelease_lib      = testrelease+'/lib/'
    testreleaseBaBar_inc = testrelease+'/BaBar/include'
else:
    testrelease      = None
    testrelease_lib  = None
    testreleaseBaBar_inc = None

# Define scons-local environment - it will be exported later.
env = Environment( CPPPATH=[ testrelease,
                             testreleaseBaBar_inc,
                             base,
                             base+'/BaBar/include',
                             framework,
                             boost_inc,
                             clhep_inc,
                             cppunit_dir+'/include',
                             heppdt_inc,
                             libsigcpp_inc+'/sigc++-2.0',
                             libsigcpp_lib+'/sigc++-2.0/include',
                             python_dir+'/include',
                             root_dir+'/include',
                             scons_dir+'/include',
                           ],
                   LIBPATH=[ testrelease_lib,
                             base+'/lib',
                             framework+'/tmp/lib',
                             boost_lib,
                             clhep_base+'/lib',
                             cppunit_dir+'/lib',
                             heppdt_lib,
                             libsigcpp_lib,
                             python_dir+'/lib',
                             root_dir+'/lib',
                             scons_dir+'/lib',
                             '/lib', '/usr/X11R6/lib',
                           ],
                   ENV={ 'PATH' : os.environ['PATH'], 
                         'LD_LIBRARY_PATH': os.environ['LD_LIBRARY_PATH'],
                       },
                   FORTRAN = 'g77',
                 )

# Define the rule for building dictionaries.
genreflex_flags = '--deep --fail_on_warnings  --capabilities=classes_ids.cc '\
                + '-DCMS_DICT_IMPL -D_REENTRANT -DGNU_SOURCE  -DGNU_GCC '\
                + '-DPROJECT_NAME="CMSSW" -DPROJECT_VERSION="CMSSW_3_0_0_pre2"'
aa="if   t1=`expr ${TARGET} : '\(.*\)_dict_plugin.cpp'`;then t2=$${t1}_map_plugin.cpp; t1=$${t1}_dict_plugin.cpp;"\
  +"elif t1=`expr ${TARGET} : '\(.*\)_map_plugin.cpp'`;then t2=$${t1}_map_plugin.cpp; t1=$${t1}_dict_plugin.cpp; fi;"\
  +"genreflex $SOURCE -s ${SOURCE.srcdir}/classes_def.xml $_CPPINCFLAGS"\
  +" -o $$t1 "\
  +genreflex_flags\
  +" && mv ${TARGET.dir}/classes_ids.cc $$t2"

genreflex = Builder(action=aa)
env.Append(BUILDERS = {'DictionarySource' : genreflex})

# Set compile and link flags.
SetOption('warn', 'no-fortran-cxx-mix')
env.MergeFlags('-g')
env.MergeFlags('-O2')
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

# Tell scons to operate on all of the SConscript files found in the previous step.
env.SConscript(ss)

# This tells emacs to view this file in python mode.
# Local Variables:
# mode:python
# End:
