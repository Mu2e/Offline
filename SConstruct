#!/usr/bin/env python
#
# Build a Mu2e base release or test release.
#
# $Id: SConstruct,v 1.54 2014/08/02 05:23:26 gandr Exp $
# $Author: gandr $
# $Date: 2014/08/02 05:23:26 $
#
# Original author Rob Kutschke.
#
import os, re, string
import sys

# Check that the release-specific setup has been run.
if not os.environ.has_key('MU2E_BASE_RELEASE'):
    sys.exit('You must setup a Mu2e base release before running scons.\nExiting.')

# Tell scons about a new command line option that controls the selection of compiler and linker flags.
AddOption('--mu2elevel',
          dest='mu2elevel',
          type='string',
          nargs=1,
          action='store',
          metavar='DIR',
          default='prof',
          help='Select debug build')

# Tell scons about a new command line option that controls the selection of graphics
AddOption('--mu2egs',
          dest='mu2egs',
          type='string',
          nargs=1,
          action='store',
          metavar='DIR',
          default='OGL',
          help='Select graphics system')

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
root_inc      = os.environ['ROOT_INC']
root_sys      = os.environ['ROOTSYS']
fhicl_inc     = os.environ['FHICLCPP_INC']
fhicl_lib     = os.environ['FHICLCPP_LIB']
sqlite_inc     = os.environ['SQLITE_INC']
sqlite_lib     = os.environ['SQLITE_LIB']
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
babarlibs = [ 'mu2e_BaBar_KalmanTrack',     'mu2e_BaBar_DetectorModel',      'mu2e_BaBar_TrkBase',    'mu2e_BaBar_BField',
              'mu2e_BaBar_TrajGeom',        'mu2e_BaBar_BbrGeom',            'mu2e_BaBar_difAlgebra', 'mu2e_BaBar_ProbTools',
              'mu2e_BaBar_BaBar',           'mu2e_BaBar_CLHEP_src_Geometry', 'mu2e_BaBar_MatEnv',
              'mu2e_BaBar_Dch_DchGeomBase', 'mu2e_BaBar_Dch_DchGeom' ]

# Define scons-local environment - it will be exported later.
osenv = {}
for var in [ 'LD_LIBRARY_PATH',  'GCC_FQ_DIR',  'PATH', 'PYTHONPATH',  'ROOTSYS', 'PYTHON_ROOT', 'PYTHON_DIR', 'SQLITE_DIR', 'SQLITE_FQ_DIR' ]:
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
                             sqlite_inc,
                             cetlib_inc,
                             cpp0x_inc,
                             boost_inc,
                             clhep_inc,
                             cppunit_dir+'/include',
                             heppdt_inc,
                             root_inc,
                             xercesc_inc,
                           ],
                   LIBPATH=[ libpath_frag,
                             base+'/lib',
                             art_lib,
                             mesfac_lib,
                             fhicl_lib,
                             sqlite_lib,
                             cetlib_lib,
                             cpp0x_lib,
                             boost_lib,
                             clhep_lib,
                             cppunit_dir+'/lib',
                             heppdt_lib,
                             root_sys+'/lib',
                             '/lib', '/usr/X11R6/lib',
                             xercesc_root+'/lib',
                           ],
                   ENV=osenv,
                   FORTRAN = 'gfortran',
                   BABARLIBS = [ babarlibs ]
                 )

# Define the rule for building dictionaries.
genreflex = Builder(action="./genreflex.sh $SOURCE $TARGET  \"$_CPPINCFLAGS\"")
env.Append(BUILDERS = {'DictionarySource' : genreflex})

# Get the flag that controls compiler options. Check that it is legal.
# There is probably a way to tell AddOption to do this test internally.
level = GetOption('mu2elevel')
known_levels = ['prof', 'debug' ]
if not level in known_levels:
    print 'Unrecognized value for --mu2elevel ' + level
    print '   The value must be one of the known levels: '  + str(known_levels)
    raise Exception('foo')

graphicssys = GetOption('mu2egs')
known_gs = ['OGL', 'Qt' ]
if not graphicssys in known_gs:
    print 'Unrecognized value for --mu2egs ' + graphicssys
    print '   The value must be one of the known systems: ' + str(known_gs)
    raise Exception('gs')

# the following may not be needed as the compiler rebuilds it all, but
# we may still want to "latch" to qt

qtfilename='.isqt'
if os.path.exists(qtfilename):
    qtf = open(qtfilename,'r')
    rgs = qtf.readline(100)
    rgs.strip()
    if graphicssys!=rgs:
        print 'Inconsitent build; the previous --mu2egs was: ' \
            + str(rgs) + ' current one is ' + graphicssys
        print 'inspect (remove?) file: ' +  qtfilename + ' or verify option --mu2egs'
        raise Exception('gs')
else:
    if graphicssys == 'Qt':
        qtf = open(qtfilename,'w')
        qtf.write(graphicssys)

env.Append( MU2EOPTS = [level, graphicssys] );

# Set compile and link flags.
SetOption('warn', 'no-fortran-cxx-mix')
env.MergeFlags('-std=c++1y')
env.MergeFlags('-rdynamic')
env.MergeFlags('-Wall')
env.MergeFlags('-Wno-unused-local-typedefs')
env.MergeFlags('-g')
if level == 'prof':
    env.MergeFlags('-O3')
    env.MergeFlags('-fno-omit-frame-pointer')
    env.MergeFlags('-DNDEBUG')

if level == 'debug':
    env.MergeFlags('-O0')

# Extract gcc version.  Some libraries have this version embedded in their names.
ff = os.popen('g++ --version'); ll = ff.readline(); ff.close()
gcc_version = ll[10:13]
env.gcc_ver=gcc_version.replace('.','')

# This comes from: root-config --cflags --glibs
# Then guess at the correct location of Spectrum and MLP.
rootlibs = [ 'Core', 'Cint', 'RIO', 'Net', 'Hist', 'Spectrum', 'MLP', 'Graf', 'Graf3d', 'Gpad', 'Tree',
             'Rint', 'Postscript', 'Matrix', 'Physics', 'MathCore', 'Thread', 'Gui', 'm', 'dl' ]
env.Append( ROOTLIBS = rootlibs );

# Make the modified environment visible to all of the SConscript files
Export('env')

# Walk the directory tree to locate all SConscript files.
ss=[]
for root,dirs,files in os.walk('.'):
    for file in files:
        if file == 'SConscript': ss.append('%s/%s'%(root[2:],file))
        pass
    pass

# Define a helper class to construct names of .so libaries. Make an instance of it available to the SConscript files.
class mu2e_helper:
    """mu2e_helper: class to produce library names"""
#   This appears to behave like c++ static member and is initialized at class defintion time.
    sourceroot =  os.path.abspath('.')
#
#   Accesor
#
    def base(self):
        return self.sourceroot
#
#   Build the name of the shared library into which non-plugin compiled code will be inserted.
#   Two versions: with and without the '#/lib' path prefix.
#
    def libname(self):
        relpath = os.path.relpath('.',self.sourceroot)
        tokens = string.split(relpath,'/')
        if len(tokens) > 1:
            if tokens[len(tokens)-1] == 'src':
                tokens.pop()
                pass
            pass
        return 'mu2e_' + string.join(tokens,'_')
    def prefixed_libname(self):
        return '#/lib/' + self.libname()
#
#   Build the name of the shared library into which plugin code will be inserted.
#   Two versions: with and without the '#/lib' path prefix.
#
    def plugin_libname(self,sourcename):
        return self.libname() + '_' + sourcename[:sourcename.find('.cc')]
    def prefixed_plugin_libname(self,sourcename):
        return '#/lib/' + self.plugin_libname(sourcename)
#
#   Build a list of plugins to be biult.
#
    def plugin_cc(self):
        return Glob('*_module.cc', strings=True) + Glob('*_service.cc', strings=True) + Glob('*_source.cc', strings=True)
#
#   Build a list of .cc files that are not plugings; these go into the library named after the directory.
#
    def non_plugin_cc(self):
        tmp = non_plugin_cc = Glob('*.cc', strings=True)
        for cc in self.plugin_cc(): tmp.remove(cc)
        return tmp
#
#   Names need to build the _dict and _map libraries.
#
    def dict_tmp_name(self):
        relpath = os.path.relpath('.',self.sourceroot)
        return '#/tmp/src/' + relpath + '/' + self.libname() + '_dict.cpp'

    def map_tmp_name(self):
        relpath = os.path.relpath('.',self.sourceroot)
        return '#/tmp/src/' + relpath + '/' + self.libname() + '_map.cpp'

    def dict_libname(self):
        relpath = os.path.relpath('.',self.sourceroot)
        return self.libname() + '_dict'

    def map_libname(self):
        relpath = os.path.relpath('.',self.sourceroot)
        return self.libname() + '_map'

    def prefixed_dict_libname(self):
        return '#/lib/' + self.dict_libname()

    def prefixed_map_libname(self):
        return '#/lib/' + self.map_libname()
#
#   Make the main library.
#
    def make_mainlib( self, userlibs, cppf=[], pf=[], addfortran=False ):
        non_plugin_cc = self.non_plugin_cc()
        if addfortran:
            fortran = Glob('*.f', strings=True)
            non_plugin_cc = [ non_plugin_cc, fortran]
            pass
        libs = []
        if non_plugin_cc:
            env.SharedLibrary( self.prefixed_libname(),
                               non_plugin_cc,
                               LIBS=[ userlibs ],
                               CPPFLAGS=cppf,
                               parse_flags=pf
                              )
            libs = [ self.libname() ]
            pass
        return libs
#
#   Make one plugin library ( but does not work for _dict and _map plugins )
#
    def make_plugin( self, cc, userlibs, cppf = [], pf = []):
        env.SharedLibrary( self.prefixed_plugin_libname(cc),
                           cc,
                           LIBS=[ userlibs, ],
                           CPPFLAGS=cppf,
                           parse_flags=pf
                           )
#
#   Make all plugin libraries, excluding _dict and _map; this works if all libraries need the same link list.
#
    def make_plugins( self, userlibs, exclude_cc = [], cppf = [], pf = [] ):
        plugin_cc = self.plugin_cc()
        for cc in exclude_cc: plugin_cc.remove(cc)
        for cc in plugin_cc:
            env.SharedLibrary( self.prefixed_plugin_libname(cc),
                               cc,
                               LIBS=[ userlibs ],
                               CPPFLAGS=cppf,
                               parse_flags=pf
                               )

#
#   Make the dictionary and map plugins.
#
    def make_dict_and_map( self, userlibs ):
        if os.path.exists('classes.h'):
            if os.path.exists('classes_def.xml'):
                env.DictionarySource([ self.dict_tmp_name(),
                                       self.map_tmp_name() ],
                                     [ 'classes.h', 'classes_def.xml'] )
                env.SharedLibrary( self.prefixed_dict_libname(),
                                   self.dict_tmp_name(),
                                   LIBS=[ userlibs ]
                                   )
                env.SharedLibrary( self.prefixed_map_libname(),
                                   self.map_tmp_name()
                                   )

# Export the class so that it can be used in the SConscript files
# For reasons I don't understand, this must come before the env.SConscript(ss) line.
Export('mu2e_helper')

# Tell scons to operate on all of the SConscript files found by walking the directory tree.
env.SConscript(ss)

# for -c remove all files from ./tmp and ./lib directories to avoid stale files
cleanopt = GetOption('clean')
if (cleanopt and not COMMAND_LINE_TARGETS):
    print "running with -c and no specific target"
    print "will also remove all files in tmp and lib directories"

    for top, dirs, files in os.walk("./lib"):
        for name in files:
            ff =  os.path.join(top, name)
            print "removing file ", ff
            os.unlink (ff)

    for top, dirs, files in os.walk("./tmp"):
        for name in files:
            ff =  os.path.join(top, name)
            print "removing file ", ff
            os.unlink (ff)

    for ff in ("EventDisplay/src/EventDisplayDict.cc", "EventDisplay/src/EventDisplayDict.h"):
        print "removing file ", ff
        os.unlink (ff)

# This tells emacs to view this file in python mode.
# Local Variables:
# mode:python
# End:
