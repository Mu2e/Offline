#!/usr/bin/env python
################################################################################
#           HOW RUN THE TRIGGER-FCL GENERATOR SCRIPT                           #
#------------------------------------------------------------------------------#
# 
# Trigger/python/genTriggerFcl.py -c Trigger/data/allTrig.config 
# or just
# Trigger/python/genTriggerFcl.py -c allTrig
# add "-o" to create online main fcl
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

def appendEpilog(trig_path, projectDir, verbose, doWrite, sourceFiles, targetFiles):

    trk_filters      = ['EventPrescale','SDCountFilter','TCFilter', 'HSFilter', 'TSFilter','Prescale']
    helix_filters    = ['EventPrescale','SDCountFilter','TCFilter', 'HSFilter', 'Prescale']
    tc_filters       = ['EventPrescale','SDCountFilter','TCFilter', 'Prescale']
    calo_filters     = ['EventPrescale','CDCountFilter','Filter'  , 'Prescale']
    unbiased_filters = ['Prescale']
    minbias_filters  = ['EventPrescale','Filter'       , 'Prescale']

    filters     = []

    #understand which kind of trigger path are we dealing with
    if "Seed" in trig_path:
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
    subEpilogDirName = projectDir + "/" + trig_path

    if verbose :
        print("Creating directory {}".format(subEpilogDirName))
    if doWrite :
        if not os.path.exists(subEpilogDirName) :
            os.makedirs(subEpilogDirName)

    subEpilogName = subEpilogDirName + ".fcl"
    if verbose :
        print("Creating {}".format(subEpilogName))
    if doWrite :
        subEpilogFile = open(subEpilogName,"w")
    targetFiles.append(subEpilogName)

    for filter in filters :
        filterName       = trig_path+filter

        subSubEpilogInputFileName = "Trigger/data/" + trig_path + "/main_"+ filterName + '.fcl' 
        sourceFiles.append(subSubEpilogInputFileName)
        subSubEpilogFileName = subEpilogDirName + "/main_"+ filterName + '.fcl' 
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

        epilog=("\n#include \""+subSubEpilogFileName +"\"")

        if doWrite :
            subEpilogFile.write(epilog)

    if doWrite:
        subEpilogFile.close()

    # return a line to be added to the main epilog file
    # so it can include the files we just wrote
    subEpilogInclude=("\n#include \""+subEpilogName+"\"")

    return subEpilogInclude

#
# routine that drives the fcl generation
#
# returns the list of files input and output, for use in scons
#

def generate(configFileText="allTrig", online=False, verbose=True, doWrite=True):

    if verbose :
        print("doWrite = {}".format(doWrite))

    # when we run from SConscript, the cwd is the python subdir
    # but all file name are relative to Offline, so go there
    cwd = os.getcwd()
    words = cwd.split("/")
    if words[-1] == "python" :
        os.chdir("../..")
        #print os.getcwd()

    # lists of files to send to scons for dependencies
    sourceFiles = []
    targetFiles = []

    # allow "Trigger/data/allTrig.config" or "allTrig.config" or "allTrig"
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
    configFileName = "Trigger/data/" + configFileBaseName + ".config"
    sourceFiles.append(configFileName)

    trig_prolog_files = [
        'Trigger/fcl/templates.fcl',
        'TrkFilters/fcl/prolog_trigger.fcl',
        'CaloFilters/fcl/prolog_trigger.fcl'
        ]
    for fn in trig_prolog_files:
        sourceFiles.append(fn)

    hasFilteroutput = False
    
    projectDir = "gen/fcl/Trigger"
    if online :
        projectDir = projectDir + "/online"
    else :
        projectDir = projectDir + "/offline"
    projectDir = projectDir + "/" + configFileBaseName

    if doWrite :
        if not os.path.exists(projectDir) :
            os.makedirs(projectDir)

    # this is the main fcl file, start with copy from template
    mainFclFileName = projectDir+"/main.fcl"
    targetFiles.append(mainFclFileName)
    if verbose :
        print("Creating {}".format(mainFclFileName))

    if online :
        templateFileName = "Trigger/fcl/main_online.fcl"
    else:
        templateFileName = "Trigger/fcl/main.fcl"
    sourceFiles.append(templateFileName)

    if doWrite :
        copyfile(templateFileName, mainFclFileName)

    # now open main fcl file for appending paths
    if doWrite :
        mainFclFile = open(mainFclFileName,"a",encoding="utf-8")

    path_list = ""
    trig_list = ""

    mainEpilogFileName   = projectDir + "/" + "allPaths.fcl"
    targetFiles.append(mainEpilogFileName)
    if verbose :
        print("Creating {}".format(mainEpilogFileName))
    if doWrite :
        mainEpilogFile   = open(mainEpilogFileName, "w");

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

            digi_path = ""
            if online :
                if 'tpr'  in pathName or 'cpr' in pathName or 'Sd' in pathName: 
                    digi_path += "makeSD, "
                if 'calo' in pathName or 'cpr' in pathName or 'Cd' in pathName: 
                    digi_path += "CaloDigiFromShower, "

            new_path = ("\nphysics."+pathName+"_trigger"+" : [ "+ digi_path +"@sequence::Trigger.paths."+pathName+" ] \n")

            #now append the epilog files for setting the filters in the path
            subEpilogInclude = appendEpilog(pathName, projectDir, verbose, 
                                            doWrite, sourceFiles, targetFiles)

            if doWrite :
                mainEpilogFile.write(subEpilogInclude)
                mainEpilogFile.write(new_path)

        else:
            # triggerOutput keyword means create an output path
            trigerOutput_line= ("\nphysics.out : [ readTriggerInfo, triggerOutput ]"+" \n\n")
            if doWrite :
                mainFclFile.write(trigerOutput_line)
            hasFilteroutput = True

    if doWrite :
        mainEpilogFile.write("\n")

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

    # include the main epilog file, which includes the others, in main fcl
    if doWrite :
        mainFclFile.write("\n#include \""+projectDir+"/allPaths.fcl\"\n")
        mainFclFile.close()

    if verbose :
        print("")
        print("main fcl: {}".format(mainFclFileName))
        print("Top level epilog file: ".format(mainEpilogFileName))
        print("")

    # now cd back to where we started
    os.chdir(cwd)

    #print len(sourceFiles),sourceFiles
    #print len(targetFiles),targetFiles
    return sourceFiles, targetFiles

#
# main, runs if started at the command line
#

if __name__ == "__main__":

    parser = ArgumentParser()
    parser.add_argument("-c", "--config-file", dest="configFileText",
                        help="file with Trigger configuration. Paths available are: unbiased, minimumbiasSdCount,largeSdCount, minimumbiasCdCount,largeCdCount, caloOnly, caloMixed, caloCosmicMuon, tprDeMSeed, tprDePSeed, cprDeMSeed, cprDePSeed, triggerOutput", metavar="FILE")
    parser.add_argument("-o", "--online", dest="online", action="store_true",
                        help="if present, use the online main fcl file template instead of offline")
    parser.add_argument("-q", "--quiet",
                        action="store_false", dest="verbose", default=True,
                        help="don't print status messages to stdout")
    
    args = parser.parse_args()
    if args.verbose :
        print("Config file name: {}".format(args.configFileText))
        print("Online flag: {}".format(str(args.online)))

    generate(args.configFileText, args.online, args.verbose, True)


    exit(0)
