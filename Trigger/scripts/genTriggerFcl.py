################################################################################
#           HOW RUN THE TRIGGER-FCL GENERATOR SCRIPT                           #
#------------------------------------------------------------------------------#
# 1) 
#
#
# 2) 

import re
import sys
import string
import os
import shutil
from shutil import copyfile

from argparse import ArgumentParser

def appendEpilog(trig_path, output_file, project_name, trig_path_counter):

    trk_filters      = ['EventPrescale','SDCountFilter','TCFilter', 'HSFilter', 'TSFilter','Prescale']
    helix_filters    = ['EventPrescale','SDCountFilter','TCFilter', 'HSFilter', 'Prescale']
    tc_filters       = ['EventPrescale','SDCountFilter','TCFilter', 'Prescale']
    calo_filters     = ['EventPrescale','CDCountFilter','Filter'  , 'Prescale']
    unbiased_filters = ['Prescale']
    minbias_filters  = ['EventPrescale','Filter'       , 'Prescale']

    filters     = []

#undestand which kind of trigger path are we dealing with
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

    trk_alg     = trig_path

#create the directory where to fill all the sub-epilog files
    project_dir = os.environ["MU2E_BASE_RELEASE"] +"/"+ project_name

    print project_dir
#check if the dir with the subconfig files exists, otherwise create it
    if os.path.exists(project_dir) == False:
        #os.makedirs(project_dir)
        project_name = 'TriggerEpilogs'
        project_dir  = os.environ["MU2E_BASE_RELEASE"] +"/TriggerEpilogs"

    #create the sub-epilog file
    sub_epilog_name = project_name + "/" + trig_path + ".fcl"
    sub_epilog_file = open(sub_epilog_name,"a")

    is_first = True

#    for i in range(0,len(trk_filters)):
#        filterName       = trk_alg+trk_filters[i]
    for i in range(0,len(filters)):
        filterName       = trk_alg+filters[i]
        if os.path.exists(project_dir) == False:
            os.makedirs(project_dir)

        project_subdir   = project_dir + "/"+trig_path
        subconfig_file   = project_subdir +"/"+filterName+ ".config"

        subconfig_exists = os.path.isfile(subconfig_file)
#if the file doesn't exist use the default
        epilog_fileName  = ""
        if subconfig_exists == False:
#            subconfig_file = os.environ["MU2E_BASE_RELEASE"]+'/Trigger/scripts/inputs/' + trig_path + "/main_"+ filterName + ".config"
#            epilog_fileName = 'Trigger/scripts/inputs/' + trig_path + "/main_"+ filterName + '.fcl' 
            epilog_fileName = project_name + "/" + trig_path + "/main_"+ filterName + '.fcl' 
            dirDefault      = 'Trigger/scripts/inputs/' + trig_path
            dirDest         = project_name + "/" + trig_path
            if is_first == True: 
                shutil.copytree(dirDefault, dirDest)
                is_first = False
        else:
            subconfig_file = open(subconfig_file,"r")

#check if the subdir for the epilog file exists, otherwise create it
            epilog_subdir      = os.environ["MU2E_BASE_RELEASE"] +"/"+ project_name + '/' + trig_path
            if os.path.exists(epilog_subdir) == False:
                os.makedirs(epilog_subdir)
                print epilog_subdir

#create the epilog file
                epilog_fileName = project_name + '/' + trig_path +'/'+ filterName+'.fcl'
                epilog_file     = open(epilog_fileName,"a")

#now read the sub-config file
                for line in subconfig_file:
                    vec_varName  = []
                    vec_varValue = []
                    for t in line.split():
                        try:
                            vec_varValue.append(float(t))
                        except ValueError:
                            vec_varName.append(t)
                            
                            if len(vec_varValue) == 1:
                                new_line   =("physics.filters."+filterName+"."+str(vec_varName[0])+" : "+str(vec_varValue[0])+ "\n")
                                epilog_file.write(new_line)
                            elif len(vec_varName) == 2:
                                new_line   =("physics.filters."+filterName+"."+str(vec_varName[0])+" : "+str(vec_varName[1])+ "\n")
                                epilog_file.write(new_line)
                                epilog_file.close()
        # set the TriggerAlg using the trig_path_counter
        epilog_file     = open(epilog_fileName,"a")
        trigAlg_line    = ("physics.filters."+filterName+".triggerPath        " + " : " + "\""+trk_alg+"_path\" \n")
        epilog_file.write(trigAlg_line)

        epilog=("\n#include \""+epilog_fileName +"\"")
        print epilog

        sub_epilog_file.write(epilog)
    
    sub_epilog=("\n#include \""+sub_epilog_name+"\"")
    output_file.write(sub_epilog)

    return True

parser = ArgumentParser()
parser.add_argument("-c", "--config-file", dest="configfilename",
                    help="file with Trigger configuration. Paths available are: unbiased, minimumbiasSdCount,largeSdCount, minimumbiasCdCount,largeCdCount, caloOnly, caloMixed, caloCosmicMuon, tprDeMSeed, tprDePSeed, cprDeMSeed, cprDePSeed, triggerOutput", metavar="FILE")
parser.add_argument("-o", "--output-file", dest="outputfilename",
                    help="name of the generated fcl file for running the Trigger executable", metavar="OUTPUT")
parser.add_argument("-t", "--template", dest="templatefilename", default= os.environ["MU2E_BASE_RELEASE"]+'/Trigger/scripts/main.fcl',
                    help="name of the fcl template file used to create the Trigger executable", metavar="TEMPLATE")
parser.add_argument("-q", "--quiet",
                    action="store_false", dest="verbose", default=True,
                    help="don't print status messages to stdout")

args = parser.parse_args()
#print args

fh = open(args.configfilename, "r")

trig_paths = [
    #unbiased trigger path
    "unbiased",
    #minimum bias filters
    "minimumbiasSdCount",
    #path for selecting events with large ammount of strawDigis
    "largeSdCount",
    #path for the calorimeter only trigger
    "caloMVACE",
    #path for calorimeter cosmic muon calibration
    "caloCalibCosmic",
    #paths for TrkPatRec downstream e- and e+
    "tprSeedDeM",
    "tprSeedDeP",
    #paths for CalPatRec downstream e- and e+
    "cprSeedDeM",
    "cprSeedDeP"
    ]

trig_prolog_files = [
    os.environ["MU2E_BASE_RELEASE"]+'/Trigger/fcl/templates.fcl',
    os.environ["MU2E_BASE_RELEASE"]+'/TrkFilters/fcl/prolog_trigger.fcl',
    os.environ["MU2E_BASE_RELEASE"]+'/CaloFilters/fcl/prolog_trigger.fcl'
    ]

#print fh.readline()
new_file   = args.outputfilename
input_file = args.templatefilename 

isOnlineMode = False
if 'Mu2eProducer' in open(input_file).read():
    isOnlineMode = True

copyfile(input_file, new_file)

new_file = open(new_file,"a")

path_list = ""
trig_list = ""

hasFilteroutput=0

tmp_name     = args.outputfilename
fname_len    = len(tmp_name)-4
project_name = tmp_name[0:fname_len]
project_dir  = os.environ["MU2E_BASE_RELEASE"] + "/" + project_name

if os.path.exists(project_dir) == False:
    project_name = "TriggerEpilogs"
    project_dir  = os.environ["MU2E_BASE_RELEASE"] +"/" + project_name
    os.makedirs(project_dir)

new_epilog   = project_name + "/" + "allPaths.fcl"
new_epilog   = open(new_epilog, "a");

#we need a counter for associating the trig bits
#with a given string in the TriggerInfoMerger module
trigger_path_counter= 0

for line in fh:
    vec_path     = []
    vec_prescale = []
    for t in line.split():
        try:
            vec_prescale.append(int(t))
            # new_line= ("physics.filters."+str(vec_path[0])+"EventPrescale"+" : "+ str(vec_prescale[0])+" \n")
            # new_file.write(new_line)
        except ValueError:
            vec_path.append(t)

    if vec_path[0] != "triggerOutput":
        #check if the name of the path is present in the prolog_trigger files
        pathCheck=False
        for i in range(0, len(trig_prolog_files)):
            if vec_path[0] in open(trig_prolog_files[i]).read():
                pathCheck  = True
        if pathCheck == False: 
            print (vec_path[0]+" NOT FOUND IN ANY PROLOG_TRIGGER.FCL FILES. PLEASE, CHECK THE INPUT FILE PROVIDED \n")
            exit()
        print pathCheck

        if path_list == "":
            path_list += vec_path[0]+"_path"
            trig_list += "\""+vec_path[0]+"\""
        else:
            path_list += ", "+vec_path[0]+"_path"
            trig_list += ", "+"\""+vec_path[0]+"\""
    
        if isOnlineMode == False:
            new_path= ("\nphysics."+vec_path[0]+"_path"+" : [ @sequence::Trigger.paths."+vec_path[0]+" ] \n") 
        else:
            digi_path=""
            if 'tpr'  in vec_path[0] or 'cpr' in vec_path[0] or 'Sd' in vec_path[0]: 
                digi_path += "makeSD, "
            if 'calo' in vec_path[0] or 'cpr' in vec_path[0] or 'Cd' in vec_path[0]: 
                digi_path += "CaloDigiFromShower, "
            new_path= ("physics."+vec_path[0]+"_path"+" : [ "+ digi_path +"@sequence::Trigger.paths."+vec_path[0]+" ] \n")

#now append the epilog files for setting the filters in the path
        
#        if "Seed" in vec_path[0]:
        appendEpilog(str(vec_path[0]), new_epilog, project_name, trigger_path_counter)

        new_epilog.write(new_path)

    else:
        trigerOutput_line= ("\nphysics.out : [ readTriggerInfo, triggerOutput ]"+" \n\n")
        new_file.write(trigerOutput_line)
        hasFilteroutput=1

    trigger_path_counter = trigger_path_counter + 1
    # if vec_path[0] != "triggerOutput":
    #     if len(vec_prescale) == 1:
    #         new_line= ("physics.filters."+str(vec_path[0])+"Prescale.nPrescale"+" : "+ str(vec_prescale[0])+" \n \n")
    #         new_epilog.write(new_line)
    #     elif len(vec_prescale) == 2:
    #         new_line= ("physics.filters."+str(vec_path[0])+"Prescale.nPrescale"+" : "+ str(vec_prescale[0])+" \n")
    #         new_epilog.write(new_line)
            
    #         new_line= ("physics.filters."+str(vec_path[0])+"EventPrescale.nPrescale"+" : "+ str(vec_prescale[1])+" \n \n")
    #         new_epilog.write(new_line)

new_epilog.write("\n")
       
analyzer_line= ("physics.analyzers.readTriggerInfo.SelectEvents : [ "+path_list+" ]"+" \n")
new_file.write(analyzer_line)
analyzer_paths_list= ("physics.analyzers.readTriggerInfo.triggerPathsList : [ "+path_list+" ]"+" \n\n")
new_file.write(analyzer_paths_list)

if hasFilteroutput == 1:
    triggerOutput_line= ("outputs.triggerOutput.SelectEvents : [ "+path_list+" ]"+" \n")
    new_file.write(triggerOutput_line)

new_epilog.close()

new_file.write("\n#include \"TriggerEpilogs/allPaths.fcl\"\n")
new_file.close()
