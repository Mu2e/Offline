import re
import sys
import string
from shutil import copyfile

from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-c", "--config-file", dest="configfilename",
                    help="file with Trigger configuration. Paths available are: unbiased, minimumbiasSdCount,largeSdCount, caloOnly, caloMixed, caloCosmicMuon, tprDeMSeed, tprDePSeed, cprDeMSeed, cprDePSeed, triggerOutput", metavar="FILE")
parser.add_argument("-o", "--output-file", dest="outputfilename",
                    help="name of the generated fcl file for running the Trigger executable", metavar="OUTPUT")
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

#print fh.readline()
new_file   = args.outputfilename
input_file = 'main.fcl'
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
        if path_list == "":
            path_list += vec_path[0]+"_path"
        else:
            path_list += ", "+vec_path[0]+"_path"
    
        new_path= ("physics."+vec_path[0]+"_path"+" : [ @sequence::paths."+vec_path[0]+" ] \n") 
        new_file.write(new_path)
    else:
        trigerOutput_line= ("physics.out : [ ReadTriggerInfo, triggerOutput ]"+" \n")
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
       
analyzer_line= ("physics.analyzers.ReadTriggerInfo.SelectEvents : [ "+path_list+" ]"+" \n")
new_file.write(analyzer_line)

if hasFilteroutput == 1:
    triggerOutput_line= ("outputs.triggerOutput.SelectEvents : [ "+path_list+" ]"+" \n")
    new_file.write(triggerOutput_line)
        
new_file.close()
