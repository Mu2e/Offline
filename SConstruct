# Script to build this release.

#
# $Id: SConstruct,v 1.8 2010/09/30 14:14:27 kutschke Exp $
# $Author: kutschke $
# $Date: 2010/09/30 14:14:27 $
#
# Original author Rob Kutschke.
#
import os
import sys

if not os.environ.has_key('FRAMEWORK_DIR'):
    sys.exit('You must setup framwork before running scons.\nExiting.')

SetOption('warn', 'no-fortran-cxx-mix')

home = os.environ['FRAMEWORK_DIR']
boost_dir = os.environ['BOOST_DIR']
boost_inc = os.environ['BOOST_INC']
clhep_dir = os.environ['CLHEP_DIR']
#cmake_dir = os.environ['CMAKE_DIR']
cppunit_dir = os.environ['CPPUNIT_DIR']
gccxml_dir = os.environ['GCCXML_DIR']
heppdt_dir = os.environ['HEPPDT_DIR']
libsigcpp_dir = os.environ['LIBSIGCPP_DIR']
python_dir = os.environ['PYTHON_DIR']
root_dir = os.environ['ROOT_DIR']
scons_dir = os.environ['SCONS_DIR']

# '#' puts the current (top-level) directory into the CPPPATH.
env = Environment( CPPPATH=[ '#',
                             '#/BaBar/include',
                             '.',
                             home,
                             boost_inc,
                             clhep_dir+'/include',
                             cppunit_dir+'/include',
                             heppdt_dir+'/include',
                             libsigcpp_dir+'/include/sigc++-2.0',
                             libsigcpp_dir+'/lib/sigc++-2.0/include',
                             python_dir+'/include',
                             root_dir+'/include',
                             scons_dir+'/include',
                           ],
                   LIBPATH=[ '#/lib',
                             home+'/tmp/lib' ,
                             boost_dir+'/lib',
                              clhep_dir+'/lib',
                             cppunit_dir+'/lib',
                             heppdt_dir+'/lib',
                             libsigcpp_dir+'/lib',
                             python_dir+'/lib',
                             root_dir+'/lib',
                             scons_dir+'/lib',
                             '/lib', '/usr/X11R6/lib',
                           ],
                   ENV={ 'PATH' : os.environ['PATH'], 
                         'LD_LIBRARY_PATH': os.environ['LD_LIBRARY_PATH'],
                         'FRAMEWORK_DIR' : os.environ['FRAMEWORK_DIR'],
                       },
                   FORTRAN = 'g77',
                 )

#env.AppendUnique(CCFLAGS=['-rpath' + home+'/tmp/lib'])

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

env.MergeFlags('-g')
env.MergeFlags('-O2')
env.MergeFlags('-rdynamic')
env.MergeFlags('-Wall')

ff = os.popen('g++ --version'); ll = ff.readline(); ff.close()
gcc_version = ll[10:13]
env.gcc_ver=gcc_version.replace('.','')

Export('env')

# Find all Package/src/SConscript files
ss=[]
for root,dirs,files in os.walk('.'):
    for file in files:
        if file == 'SConscript': ss.append('%s/%s'%(root[2:],file))
        pass
    pass

# scons to bulid all of the SConscript files
env.SConscript(ss)

# This tells emacs to view this file in python mode.
# Local Variables:
# mode:python
# End:
