#!/usr/bin/env python
################################################################################
#           HOW RUN THE TRIGGER-FCL GENERATOR SCRIPT                           #
#------------------------------------------------------------------------------#
# 
# Trigger/python/genTriggerFcl.py -c Trigger/data/allPaths.config 
# or just
# Trigger/python/genTriggerFcl.py -c allPaths
#

import re
import sys
import string
import os
import shutil
from shutil import copyfile

from argparse import ArgumentParser

from codecs import open

#
# process one subdirectory (one path)
#

def appendEpilog(trig_path, relProjectDir, outDir, srcDir, verbose, doWrite, sourceFiles, targetFiles):

    trk_filters      = ['EventPrescale','SDCountFilter','TCFilter', 'HSFilter', 'TSFilter']
    helix_filters    = ['EventPrescale','SDCountFilter','TCFilter', 'HSFilter']
    tc_filters       = ['EventPrescale','SDCountFilter','TCFilter']
    calo_filters     = ['EventPrescale','CDCountFilter','Filter'  ]
    unbiased_filters = ['EventPrescale']
    minbias_filters  = ['EventPrescale','Filter'       ]
    cst_filters      = ['EventPrescale','SDCountFilter','TCFilter', 'TSFilter']
    
    filters     = []

    #understand which kind of trigger path are we dealing with
    if "Seed" in trig_path:
        if "cst" in trig_path:
            filters = cst_filters
        else:
            filters = trk_filters
    elif "Helix" in trig_path:
        filters = helix_filters
    elif "TimeCluster" in trig_path:
        filters = tc_filters
    elif "calo" in trig_path:
        filters = calo_filters
    elif "unbiased" in trig_path:
        filters = unbiased_filters
    elif "Count" in trig_path:
        filters = minbias_filters

    if len(filters) == 0 :
        print("ERROR: path {} has no associated filters".format(trig_path))
        exit(1)

    #create the sub-epilog file
    subEpilogDirName = outDir+relProjectDir + "/" + trig_path
    relSubEpilogDirName = relProjectDir + "/" + trig_path

    if verbose :
        print("Creating directory {}".format(subEpilogDirName))
    if doWrite :
        if not os.path.exists(subEpilogDirName) :
            os.makedirs(subEpilogDirName)

    subEpilogName = subEpilogDirName + ".fcl"
    relSubEpilogName = relSubEpilogDirName + ".fcl"
    if verbose :
        print("Creating {}".format(subEpilogName))
    if doWrite :
        subEpilogFile = open(subEpilogName,"w")
    targetFiles.append(subEpilogName)

    for filter in filters :
        filterName       = trig_path+filter

        subSubEpilogInputFileName = srcDir+"Trigger/data/" + trig_path + "/main_"+ filterName + '.fcl' 
        sourceFiles.append(subSubEpilogInputFileName)
        subSubEpilogFileName      = subEpilogDirName + "/main_"+ filterName + '.fcl' 
        relSubSubEpilogFileName   = relSubEpilogDirName + "/main_"+ filterName + '.fcl' 
        targetFiles.append(subSubEpilogFileName)
        if verbose:
            print("Creating {}".format(subSubEpilogFileName))

        # first copy the fcl file
        if doWrite :
            shutil.copy(subSubEpilogInputFileName,subEpilogDirName)

        # then open it and append one line
        if  doWrite :
            subSubEpilogFile = open(subSubEpilogFileName,"a")
        trigAlgLine    = ("\nphysics.filters."+filterName+".triggerPath        " + " : " + "\""+trig_path+"_trigger\" \n")
        if doWrite :
            subSubEpilogFile.write(trigAlgLine)
            subSubEpilogFile.close()

        epilog=("\n#include \""+relSubSubEpilogFileName +"\"")

        if doWrite :
            subEpilogFile.write(epilog)

    #now create the instance for the TriggerInfo Merger
    if  doWrite :
        trigInfoMergerName         = trig_path + "TriggerInfoMerger"
        subSubEpilogMergerFileName = subEpilogDirName + "/main_" + trigInfoMergerName + '.fcl'
        if verbose:
            print("Creating {}".format(subSubEpilogMergerFileName))
        subSubEpilogMergerFile     = open(subSubEpilogMergerFileName,"w+")
        subSubEpilogMergerFile.write("physics.producers."+trigInfoMergerName+" : { module_type : MergeTriggerInfo }");
        subSubEpilogMergerFile.close();

        relSubSubEpilogFileName    = relSubEpilogDirName + "/main_"+ trigInfoMergerName + '.fcl' 
        epilog=("\n#include \""+relSubSubEpilogFileName +"\"")

        subEpilogFile.write(epilog)
        subEpilogFile.close()

    # return a line to be added to the main epilog file
    # so it can include the files we just wrote
    subEpilogInclude=("\n#include \""+relSubEpilogName+"\"")

    return subEpilogInclude

#
# routine that drives the fcl generation
#
# returns the list of files input and output, for use in scons
#

def generate(configFileText="allPaths", verbose=True, doWrite=True):

    if verbose :
        print("configFileText = ",configFileText)
        print("doWrite = ",doWrite)

    # when we run from SConscript to create the targets, 
    # we are in the subdir.  When running in Offline scons, we are in Offline,
    # when in Muse, we are in Offline parent dir
    # This code makes this always run in the scons default dir, either 
    # Offline or its parent dir
    srcDir=""
    outDir=""
    owd = os.getcwd()
    if verbose:
        print ("start owd= ",owd)
    if 'MUSE_WORK_DIR' in os.environ:
        # then we are running in Muse, make adjustments
        if verbose:
            print ("running in Muse mode ")
        os.chdir(os.environ['MUSE_WORK_DIR'])
        srcDir = "Offline/"
        # like "build/sl7-prof-e20/Offline/"
        outDir = os.environ['MUSE_BUILD_BASE']+"/Offline/"
    else:
        # when we run from SConscript, the owd is the python subdir
        # but all file name are relative to Offline, so go there
        words = owd.split("/")
        if words[-1] == "python" :
            os.chdir("../..")
        # accept empty default srcDir and outDir, so build in Offline

    if verbose:
        print ("owd = ",owd)
        print ("pwd = ",os.getcwd())
        print ("srcDir = ",srcDir)
        print ("outDir = ",outDir)

    # lists of files to send to scons for dependencies
    sourceFiles = []
    targetFiles = []

    # allow "Trigger/data/allPaths.config" or "allPaths.config" or "allPaths"
    tempArr = configFileText.split("/")
    if len(tempArr) > 1:
        temp = "/".join (tempArr[:-1])
        if temp != "Trigger/data" :
            print("ERROR config file must be in directory Trigger/data")
            exit(1)
    tempArr = tempArr[-1].split(".")
    if len(tempArr) > 1:
        if tempArr[1] != "config" or len(tempArr)>2 :
            print("ERROR config file must of type .config")
            exit(1)
    configFileBaseName = tempArr[0];
    configFileName = srcDir+"Trigger/data/" + configFileBaseName + ".config"
    sourceFiles.append(configFileName)

    trig_prolog_files = [
        srcDir+'Trigger/fcl/prolog_trigger.fcl',
        srcDir+'TrkFilters/fcl/prolog_trigger.fcl',
        srcDir+'CaloFilters/fcl/prolog_trigger.fcl',
        srcDir+'CosmicReco/fcl/prolog_trigger.fcl'
        ]
    for fn in trig_prolog_files:
        sourceFiles.append(fn)

    hasFilteroutput = False
    
    relProjectDir = "gen/fcl/Trigger/" + configFileBaseName
    projectDir = outDir+relProjectDir

    if doWrite :
        if not os.path.exists(projectDir) :
            os.makedirs(projectDir)

    # this is the main fcl file, start with copy from template
    mainFclFileName = projectDir+"/main.fcl"
    targetFiles.append(mainFclFileName)
    if verbose :
        print("Creating {}".format(mainFclFileName))

    templateFileName = srcDir+"Trigger/fcl/main.fcl"
    sourceFiles.append(templateFileName)

    if doWrite :
        copyfile(templateFileName, mainFclFileName)

    # now open main fcl file for appending paths
    if doWrite :
        mainFclFile = open(mainFclFileName,"a",encoding="utf-8")

    path_list = ""
    trig_list = ""

    mainEpilogFileName   = projectDir + "/" + "{}.fcl".format(configFileBaseName)
    mainEpilogTimingFileName = projectDir + "/" + "{}_timing.fcl".format(configFileBaseName)

    targetFiles.append(mainEpilogFileName)
    if verbose :
        print("Creating {}".format(mainEpilogFileName))
    if doWrite :
        mainEpilogFile       = open(mainEpilogFileName, "w");
        mainEpilogTimingFile = open(mainEpilogTimingFileName, "w");

    #
    # main loop over lines in the config file
    #

    configFile = open(configFileName, "r")

    for line in configFile:

        line = line.strip() # strip whitespace
        if len(line)==0:
            continue  # skip empty lines

        # parse line: path [prescale] [prescale]
        words = line.split()
        pathName = words[0]

        if pathName != "triggerOutput":

            # check if the name of the path is present in the prolog_trigger files
            pathCheck=False
            for i in range(0, len(trig_prolog_files)):
                if pathName in open(trig_prolog_files[i]).read():
                    pathCheck  = True
            if pathCheck == False: 
                print ("{} NOT FOUND IN ANY PROLOG_TRIGGER.FCL FILES. PLEASE, CHECK THE INPUT FILE PROVIDED".format(pathName))
                exit(1)

            if path_list != "":
                path_list += ", "
                trig_list += ", "
            path_list += pathName+"_trigger"
            trig_list += "\""+pathName+"\""

            digi_path = "@sequence::Trigger.PrepareDigis, "
            
            new_path = ("\nphysics."+pathName+"_trigger"+" : [ "+ digi_path +"@sequence::Trigger.paths."+pathName+" ] \n")
            timing_paths = []
            if "Seed" in pathName:
                nFilters = 3
                if "cst" in pathName:
                    nFilters = 2                    
                for ind in range(nFilters):
                    timing_label = "Timing{:d}".format(ind)
                    timing_paths.append("\nphysics."+pathName+timing_label+"_trigger"+" : [ "+ digi_path +"@sequence::Trigger.paths."+pathName+timing_label+" ] \n")

            #now append the epilog files for setting the filters in the path
            subEpilogInclude = appendEpilog(pathName, relProjectDir, 
                                            outDir, srcDir, verbose, 
                                            doWrite, sourceFiles, targetFiles)

            if doWrite :
                mainEpilogFile.write(subEpilogInclude)
                mainEpilogFile.write(new_path)
                #
                mainEpilogTimingFile.write(subEpilogInclude)
                mainEpilogTimingFile.write(new_path)
                for l in range(len(timing_paths)):
                    mainEpilogTimingFile.write(timing_paths[l])
                

        else:
            # triggerOutput keyword means create an output path
            trigerOutput_line= ("\nphysics.out : [ readTriggerInfo, triggerOutput ]"+" \n\n")
            if doWrite :
                mainFclFile.write(trigerOutput_line)
            hasFilteroutput = True
    
    if doWrite :
        mainEpilogFile.write("\n")
        mainEpilogTimingFile.write("\n")

    analyzer_line= ("physics.analyzers.readTriggerInfo.SelectEvents : [ "+path_list+" ]"+" \n")
    if doWrite :
        mainFclFile.write(analyzer_line)
    analyzer_paths_list= ("physics.analyzers.readTriggerInfo.triggerPathsList : [ "+path_list+" ]"+" \n\n")
    if doWrite :
        mainFclFile.write(analyzer_paths_list)

    if hasFilteroutput :
        triggerOutput_line= ("outputs.triggerOutput.SelectEvents : [ "+path_list+" ]"+" \n")
        if doWrite :
            mainFclFile.write(triggerOutput_line)

    if doWrite :
        mainEpilogFile.close()
        mainEpilogTimingFile.close()

    # include the main epilog file, which includes the others, in main fcl
    if doWrite :
        mainFclFile.write("\n#include \""+relProjectDir+"/{}.fcl\"\n".format(configFileBaseName))
        mainFclFile.close()

    if verbose :
        print("")
        print("main fcl: {}".format(mainFclFileName))
        print("Top level epilog file: ".format(mainEpilogFileName))
        print("")

    # now cd back to where we started
    os.chdir(owd)

    return sourceFiles, targetFiles,srcDir+"Trigger/python/genTriggerFcl.py"

#
# main, runs if started at the command line
#

if __name__ == "__main__":

    parser = ArgumentParser()
    parser.add_argument("-c", "--config-file", dest="configFileText",
                        help="file with Trigger configuration. Paths available are: unbiased, minimumbiasSdCount,largeSdCount, minimumbiasCdCount,largeCdCount, caloOnly, caloMixed, caloCosmicMuon, tprDeMSeed, tprDePSeed, cprDeMSeed, cprDePSeed, triggerOutput", metavar="FILE")
    parser.add_argument("-q", "--quiet",
                        action="store_false", dest="verbose", default=True,
                        help="don't print status messages to stdout")
    
    args = parser.parse_args()
    if args.verbose :
        print("Config file name: {}".format(args.configFileText))

    generate(args.configFileText, args.verbose, True)


    exit(0)
