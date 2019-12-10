#
# This class defines procedures that most SConscript files will use
# for building. In an SConscript it is created:
#
# Import('env')
# Import('mu2e_helper')
# helper=mu2e_helper(env)
#
# Then the typical functions used are
#
# helper.make_mainlib(...)
# helper.make_plugins(...)
# helper.make_dict_and_map(...)
#

import os
import string
from glob import glob


class mu2e_helper:
    """mu2e_helper: class to produce libraries"""

    def __init__(self,env):
        self.env = env
        self.base = env['MU2EBASE']  # the top of the release
        self.codeDir = os.path.abspath('.')  # this subdir with the SConscript
        self.relpath = os.path.relpath(self.codeDir,self.base) # the difference
        # where dictionaries go: tmp/src/package/subdir
        self.tmpdir = "tmp/src/"+self.relpath
        # change string package/subdir/src to package_subdir
        tokens = self.relpath.split('/')
        if len(tokens) > 1:
            if tokens[-1] == 'src': tokens.pop()
        self.libstub = '_'.join(tokens)

        # A few places we use ClassDef in order to enable a class
        # to be fully capable at the root prompt
        # Using ClassDef forces the dictionary to be linked with the main
        # class code.  Set this True to make this happen
        self.classdef = False

    # set true if ClassDef is used, to force dictionary
    # to be linked in mainlib
    def classDef(self, tf=True):
        self.classdef = tf

    def lib_link_name(self):
        return "mu2e_"+self.libstub
    def lib_file(self):
        return "lib/libmu2e_"+self.libstub+".so"
    def plugin_lib_file(self,sourcename):
        stub = sourcename[:sourcename.find('.cc')] # file name minus the .cc
        return "lib/libmu2e_"+self.libstub + '_' + stub +".so"
    def dict_file(self):
        return self.tmpdir+"/mu2e_"+self.libstub + '_dict.cpp'
    def dict_lib_file(self):
        if self.classdef : # dictionary is in the main lib
            return "lib/libmu2e_"+self.libstub + '.so'
        else :  # dictionary is in its own lib
            return "lib/libmu2e_"+self.libstub + '_dict.so'
    def rootmap_file(self):
        return "lib/libmu2e_"+self.libstub + "_dict.rootmap"
    def pcm_file(self):
        if self.classdef : # dictionary is in the main lib
            return "lib/libmu2e_"+self.libstub + "_rdict.pcm"
        else :  # dictionary is in its own lib
            return "lib/libmu2e_"+self.libstub + "_dict_rdict.pcm"

    #
    #   Build a list of plugins to be built.
    #
    def plugin_cc(self):
        return self.env.Glob('*_module.cc', strings=True) \
            + self.env.Glob('*_service.cc', strings=True) \
            + self.env.Glob('*_source.cc', strings=True)  \
            + self.env.Glob('*_utils.cc',strings=True)    \
            + self.env.Glob('*_tool.cc',strings=True)

    #
    #   Build a list of bin source files
    #
    def bin_cc(self):
        return self.env.Glob('*_main.cc', strings=True)

    #
    #   Build a list of .cc files that are not plugings or bins;
    #   these go into the library named after the directory.
    #
    def mainlib_cc(self):
        cclist = self.env.Glob('*.cc', strings=True)
        for cc in self.plugin_cc(): cclist.remove(cc)
        for cc in self.bin_cc(): cclist.remove(cc)
        return cclist

    #
    #   Make the main library.
    #
    def make_mainlib( self, userlibs, cppf=[], pf=[], addfortran=False ):
        mainlib_cc = self.mainlib_cc()
        if addfortran:
            fortran = self.env.Glob('*.f', strings=True)
            mainlib_cc = [ mainlib_cc, fortran ]
        # if classdef is used, force dictionary into mainlib
        if self.classdef :
            mainlib_cc.append("#/"+self.dict_file())
        if mainlib_cc:
            self.env.SharedLibrary( "#/"+self.lib_file(),
                               mainlib_cc,
                               LIBS=[ userlibs],
                               CPPFLAGS=cppf,
                               parse_flags=pf
                              )
            return self.lib_link_name()
        else:
            return ""

    #
    #   Make one plugin library
    #
    def make_plugin( self, cc, userlibs, cppf = [], pf = []):
        self.env.SharedLibrary( "#/"+self.plugin_lib_file(cc),
                           cc,
                           LIBS=[ userlibs],
                           CPPFLAGS=cppf,
                           parse_flags=pf
                           )

    #
    #   Make all plugin libraries, excluding _dict and _map; this works if
    #   all libraries need the same link list.
    #
    def make_plugins( self, userlibs, exclude_cc = [], cppf = [], pf = [] ):
        plugin_cc = self.plugin_cc()
        for cc in exclude_cc: plugin_cc.remove(cc)
        for cc in plugin_cc:
            self.make_plugin(cc,userlibs, cppf, pf)

    #
    #   Make the dictionary and rootmap plugins.
    #
    def make_dict_and_map( self, userlibs=[], pf_dict=[] ):
        if os.path.exists('classes.h') and os.path.exists('classes_def.xml'):
            sources = ["classes.h","classes_def.xml"]
            targets = ["#/"+self.dict_file(),
                       "#/"+self.rootmap_file(),
                       "#/"+self.pcm_file() ]
            dflag = ""
            if self.env["BUILD"] == "debug":
                dflag = ""
            else:
                dflag = "-DNDEBUG"
            self.env.DictionarySource( targets, sources ,
                                       LIBTEXT=self.dict_lib_file(),
                                       DEBUG_FLAG=dflag)
            # if classdef is used, do not make the dictionary into its own lib,
            # it will be put in the mainlib
            if self.classdef :
                return
            # make lib for the dictionary
            self.env.SharedLibrary( "#/"+self.dict_lib_file(),
                                    "#/"+self.dict_file(),
                                    LIBS=[ userlibs ],
                                    parse_flags=pf_dict
                                    )
    #
    #   Make a bin based on binname_main.cc -> binname
    #
    def make_bin( self, target, userlibs=[], otherSource=[], bindir="" ):
        if not bindir : bindir = self.env['BINDIR']
        sourceFiles = [ target+"_main.cc" ] + otherSource
        self.env.Program(
            target = bindir+"/"+target,
            source = sourceFiles,
            LIBS   = userlibs
            )

    #
    #   Make any combo of files
    #
    def make_generic( self, source, target, command ):
        topSources = []
        topTargets = []
        for s in source :
            topSources.append("#/"+s)
        for t in target :
            topTargets.append("#/"+t)
        self.env.GenericBuild( topTargets, topSources, COMMAND=command)
