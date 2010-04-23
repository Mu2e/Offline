# Script to build this release.

#
# $Id: SConstruct,v 1.3 2010/04/23 04:12:09 kutschke Exp $
# $Author: kutschke $
# $Date: 2010/04/23 04:12:09 $
#
# Original author Rob Kutschke.
#
import os
import sys

if not os.environ.has_key('MU2E_EXTERNALS'):
    sys.exit('You must define MU2E_EXTERNALS before running scons.\nExiting.')

home = os.environ['MU2E_HOME']
externals = os.environ['MU2E_EXTERNALS']
rootdir = os.environ['ROOT_DIR']
heppdtdir = os.environ['HEPPDT_DIR']

# '#' puts the current (top-level) directory into the CPPPATH.
env = Environment( CPPPATH=[ '#',
                             '#/BaBar/include',
			     '.',
			     home,
			     externals+'/include',
                             heppdtdir+'/include',
			     rootdir+'/include/root',
			     externals+'/include/sigc++-2.0',
			     externals+'/lib/sigc++-2.0/include',
			   ],
		   LIBPATH=[ '#/lib',
			     home+'/tmp/lib' ,
			     externals+'/lib',
                             heppdtdir+'/lib',
			     rootdir+'/lib/root',
                             '/lib', '/usr/X11R6/lib',
                           ],
		   ENV={ 'PATH' : os.environ['PATH'], 
			 'LD_LIBRARY_PATH': os.environ['LD_LIBRARY_PATH'],
			 'MU2E_HOME' : os.environ['MU2E_HOME'],
		       }
                 )

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
env.MergeFlags('-rdynamic')

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
