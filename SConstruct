#
# Top driver script for scons to build a Mu2e base release or satellite release.
#
import os, re, string, sys

# Functions that do small tasks and build lists
import sconstruct_helper as sch
# handles how the input files are collected and the output files are named
from mu2e_helper import mu2e_helper

# this will contain global config info about the build
mu2eOpts = {}

# add a mu2e debug print option like "--mu2ePrint=5"
AddOption('--mu2ePrint', dest='mu2ePrint', 
          type='int',nargs=1,default=1,
          help='mu2e print level (0-10) default=1')
mu2ePrint = GetOption("mu2ePrint")
mu2eOpts["mu2ePrint"] = mu2ePrint

# add an option to print only short lines for each target
AddOption('--mu2eCompactPrint', dest='mu2eCompactPrint', 
          action="store_true",default=False,
          help='print only a short text line for each target')
mu2eCompactPrint = GetOption("mu2eCompactPrint")

mu2eOpts["mu2eCompactPrint"] = mu2eCompactPrint

# Check that some important environment variables have been set; 
# result is a dictionary of the options
moreOpts = sch.mu2eEnvironment()
mu2eOpts.update(moreOpts)

if mu2ePrint > 1:
    print ("mu2eOpts:")
    print (mu2eOpts)


if mu2ePrint > 5:
    print ("building Evironment:")
    print ("\nCPPPATH = ",sch.cppPath(mu2eOpts))
    print ("\nLIBPATH = ",sch.libPath(mu2eOpts))
    print ("\nENV = ",sch.exportedOSEnvironment())
    print ("\nFORTRAN = 'gfortran'")
    print ("\nBABARLIBS = ", sch.BaBarLibs())
    print ("\nmerge Flags =",sch.mergeFlags(mu2eOpts))

# if requested, simplfy the text printed per target
cccomstr = ""
linkcomstr = ""
genreflexcomstr = ""
if mu2eCompactPrint :
    cccomstr = "Compiling $SOURCE"
    linkcomstr = "Linking $TARGET"
    genreflexcomstr = "genreflex ${SOURCES[1]}"

# this the scons object which contains the methods to build code
env = Environment( CPPPATH = sch.cppPath(mu2eOpts),   # $ART_INC ...
                   LIBPATH = sch.libPath(mu2eOpts),   # /lib, $ART_LIB ...
                   ENV = sch.exportedOSEnvironment(), # LD_LIBRARY_PATH, ROOTSYS, ...
                   FORTRAN = 'gfortran',
                   BABARLIBS = sch.BaBarLibs(),
                   CXXCOMSTR = cccomstr,
                   SHCXXCOMSTR = cccomstr,
                   LINKCOMSTR = linkcomstr,
                   SHLINKCOMSTR= linkcomstr,
)

# Make the Compilation DB generator available in the environment
env.Tool('compilation_db')

# Define and register the rule for building dictionaries.
# sources are classes.h, classes_def.xml, 
# targets are dict.cpp, .rootmap and .pcm
# LIBTEXT is the library for the dict - not a target, only text for names
genreflex = Builder(action=Action("export HOME="+os.environ["HOME"]+"; "+"genreflex ${SOURCES[0]} -s ${SOURCES[1]} $_CPPINCFLAGS -l $LIBTEXT -o ${TARGETS[0]} --fail_on_warnings --rootmap-lib=$LIBTEXT  --rootmap=${TARGETS[1]} $DEBUG_FLAG",genreflexcomstr))
env.Append(BUILDERS = {'DictionarySource' : genreflex})

# a generic builder, some files transform to others
generic = Builder(action="$COMMAND" )
env.Append(BUILDERS = {'GenericBuild' : generic})


# this sets the build flags, like -std=c++14 -Wall -O3, etc
SetOption('warn', 'no-fortran-cxx-mix')
env.MergeFlags( sch.mergeFlags(mu2eOpts) )
# env.MergeFlags( '-lpthread' )

# env construction variables, in SConscript: var=env['VARNAME']
env.Append( ROOTLIBS = sch.rootLibs() )
env.Append( BABARLIBS = sch.BaBarLibs() )
env.Append( MU2EBASE = mu2eOpts["base"] )
env.Append( BINDIR = mu2eOpts['bindir'] )
env.Append( BUILD = mu2eOpts["build"] )
env.Append( G4VIS = mu2eOpts["g4vis"] )
env.Append( G4MT = mu2eOpts["g4mt"] )

# make the scons environment visible to all SConscript files (Import('env'))
Export('env')

# Export the class so that it can be used in the SConscript files
# comes before the env.SConscript(ss) line so it is available to SConscripts
# it is undocumented how you can export a class name, not a variable
Export('mu2e_helper')

# the list of SConscript files in the directory tree
ss = sch.sconscriptList(mu2eOpts)

# make sure lib, bin and tmp are there
sch.makeSubDirs(mu2eOpts)

# Allow a compile_commands.json to be generated if requested (scons -Q compiledb).
compileCommands = env.CompilationDatabase('compile_commands.json')
compileDb = env.Alias("compiledb", compileCommands)

# operate on the SConscript files
# regular python commands like os.path() are executed immediately as they are encontered, 
# scons builder commands like env.SharedLibrary are examined for dependences and scheduled
# to be executed in parallel, as possible
env.SConscript(ss)

#  with -c, scons will remove all dependant files it knows about.
#  this code removes orphan files caused by a parent that was removed
if ( GetOption('clean') and not COMMAND_LINE_TARGETS):
    sch.extraCleanup()

#
# create targets that are only built on demand
#

# make a copy of env with the full build environmentals so it can run exe's
build_env = env.Clone()
build_env.Append(ENV=os.environ)
# a text file of how packages depend on each other: run with 'scons DEPS'
pts = []
pts = pts + sch.PhonyTarget(build_env,'DEPS','gen/txt/deps.txt',
                            'scripts/build/bin/procs.sh DEPS')
# make the gdml file
pts = pts + sch.PhonyTarget(build_env,'GDML','gen/gdml/mu2e.gdml',
                            'scripts/build/bin/procs.sh GDML' )
# run root fast overlaps check
pts = pts + sch.PhonyTarget(build_env,'ROVERLAPS',[],
                            'scripts/build/bin/procs.sh ROVERLAPS' )
Depends(pts[-1], pts[-2]) # root check requires gdml
# run g4test_03
pts = pts + sch.PhonyTarget(build_env,'TEST03',[],
                            'scripts/build/bin/procs.sh TEST03' )
# pack git records
pts = pts + sch.PhonyTarget(build_env,'GITPACK',[],
                            'scripts/build/bin/procs.sh GITPACK' )
# remove tmp and intermediate .so files
pts = pts + sch.PhonyTarget(build_env,'RMSO',[],
                            'scripts/build/bin/procs.sh RMSO' )
# a standard validation file (~20min to run)
pts = pts + sch.PhonyTarget(build_env,'VAL0','gen/val/ceSimReco_5000.root',
                'scripts/build/bin/procs.sh VAL0' )
# combine the targets into one
pt_rel = env.AlwaysBuild(env.Alias('RELEASE', [], '#echo completed release targets'))
Depends(pt_rel, pts)


# This tells emacs to view this file in python mode.
# Local Variables:
# mode:python
# End:
# vi:syntax=python
