import re
import sys
import string
import os
from shutil import copyfile

from argparse import ArgumentParser

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
    "caloOnly",
    #path for calorimeter cosmic muon calibration
    "caloCosmicMuon",
    #paths for TrkPatRec downstream e- and e+
    "tprDeMSeed",
    "tprDePSeed",
    #paths for CalPatRec downstream e- and e+
    "cprDeMSeed",
    "cprDePSeed"
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

hasFilteroutput=0

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
        else:
            path_list += ", "+vec_path[0]+"_path"
    
        if isOnlineMode == False:
            new_path= ("physics."+vec_path[0]+"_path"+" : [ @sequence::paths."+vec_path[0]+" ] \n") 
        else:
            digi_path=""
            if 'tpr'  in vec_path[0] or 'cpr' in vec_path[0] or 'Sd' in vec_path[0]: 
                digi_path += "makeSD, "
            if 'calo' in vec_path[0] or 'cpr' in vec_path[0] or 'Cd' in vec_path[0]: 
                digi_path += "CaloDigiFromShower, "
            new_path= ("physics."+vec_path[0]+"_path"+" : [ "+ digi_path +"@sequence::paths."+vec_path[0]+" ] \n")

        new_file.write(new_path)
    else:
        trigerOutput_line= ("physics.out : [ readTriggerInfo, triggerOutput ]"+" \n")
        new_file.write(trigerOutput_line)
        hasFilteroutput=1


    if vec_path[0] != "triggerOutput":
        if len(vec_prescale) == 1:
            new_line= ("physics.filters."+str(vec_path[0])+"Prescale.nPrescale"+" : "+ str(vec_prescale[0])+" \n \n")
            new_file.write(new_line)
        elif len(vec_prescale) == 2:
            new_line= ("physics.filters."+str(vec_path[0])+"Prescale.nPrescale"+" : "+ str(vec_prescale[0])+" \n")
            new_file.write(new_line)
            
            new_line= ("physics.filters."+str(vec_path[0])+"EventPrescale.nPrescale"+" : "+ str(vec_prescale[1])+" \n \n")
            new_file.write(new_line)

new_file.write("\n")
       
analyzer_line= ("physics.analyzers.readTriggerInfo.SelectEvents : [ "+path_list+" ]"+" \n")
new_file.write(analyzer_line)

if hasFilteroutput == 1:
    triggerOutput_line= ("outputs.triggerOutput.SelectEvents : [ "+path_list+" ]"+" \n")
    new_file.write(triggerOutput_line)
        
new_file.close()
