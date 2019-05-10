from string import Template
import sys
import random
import os
import glob
import ROOT
from normalizations import *
import subprocess

dirname = sys.argv[1]
fin = open(dirname + "/settings")
lines = fin.readlines()

livetime = float(lines[0])
rue = float(lines[1])
rup = float(lines[2])
dem_emin = float(lines[3])
dep_emin = float(lines[4])
tmin = float(lines[5])
kmax = float(lines[6])

line = raw_input("Cosmic dem starting event (run:subrun:event)? ")
dem_cosmic_run = int(line.split(":")[0])
dem_cosmic_subrun = int(line.split(":")[1])
dem_cosmic_event = int(line.split(":")[2])
line = raw_input("Cosmic dep starting event (run:subrun:event)? ")
dep_cosmic_run = int(line.split(":")[0])
dep_cosmic_subrun = int(line.split(":")[1])
dep_cosmic_event = int(line.split(":")[2])
use_dayabay = 0
if len(sys.argv) > 2:
  use_dayabay = int(sys.argv[2])

ROOT.gRandom.SetSeed(0)

norms = {
  "dem-DIOLeadingLogMix": dio_normalization(livetime,dem_emin),
  "dep-DIOLeadingLogMix": 0,
  "dem-CeMLeadingLogMix": ce_normalization(livetime,rue),
  "dep-CeMLeadingLogMix": 0,
  "dem-CePLeadingLogMix": 0, 
  "dep-CePLeadingLogMix": ce_normalization(livetime,rup),
  "dem-RMCexternalMix": rmc_normalization(livetime, dem_emin, kmax, False),
  "dep-RMCexternalMix": rmc_normalization(livetime, dep_emin, kmax, False),
  "dem-RMCinternalMix": rmc_normalization(livetime, dem_emin, kmax, True),
  "dep-RMCinternalMix": rmc_normalization(livetime, dep_emin, kmax, True),
  "dem-RPCexternalMix": 1.39244803e+08*livetime, 
  "dep-RPCexternalMix": 1.59222825e+08*livetime,
  "dem-RPCinternalMix": 9.60742e+05*livetime,
  "dep-RPCinternalMix": 1.098685+06*livetime
  }

if use_dayabay:
  norms["dem-DYBCosmicMix"] = 5494.24*24*60.*60.*365*livetime
  norms["dep-DYBCosmicMix"] = 5494.24*24*60.*60.*365*livetime
else:
  norms["dem-CRYCosmicMix"] = 253400*24*60.*60.*365*livetime
  norms["dep-CRYCosmicMix"] = 253400*24*60.*60.*365*livetime

mean_reco_events = {}
filenames = {}

for signal in norms: 
  filenames[signal] = raw_input(signal + "filename? ")
  fin = ROOT.TFile(filenames[signal])
  te = fin.Get("Events")
  reco_events = te.GetEntries()
  t = fin.Get("SubRuns")
  h = ROOT.TH1F("h","h",1,0,1000000)
  t.Project("h","mu2e::GenEventCount_genCounter__ensembleMix.product()->count()")
  gen_events = h.GetMean()
  
  print "TOTAL:",gen_events,reco_events,"Efficiency:",reco_events/float(gen_events)
  mean_gen_events = norms[signal]
  mean_reco_events[signal] = mean_gen_events*reco_events/float(gen_events)
  #sample_reco_events[signal] = ROOT.gRandom.Poisson(mean_reco_events)
  #print "# events in sample:",sample_reco_events,"out of",mean_reco_events,"expected"
  print "Expect",mean_reco_events


total_sample_events = ROOT.gRandom.Poisson(sum(mean_reco_events.values()))
print "Poisson:",total_sample_events,sum(mean_reco_events.values())

fin = open("JobConfig/ensemble/SamplingInput.fcl")
t = Template(fin.read())
event_num = 0
starting_event_num = {signal: [0,0,0] for signal in norms}
if use_dayabay:
  starting_event_num["dem-DYBCosmicMix"] = [dem_cosmic_run,dem_cosmic_subrun,dem_cosmic_event]
  starting_event_num["dep-DYBCosmicMix"] = [dep_cosmic_run,dep_cosmic_subrun,dep_cosmic_event]
else:
  starting_event_num["dem-CRYCosmicMix"] = [dem_cosmic_run,dem_cosmic_subrun,dem_cosmic_event]
  starting_event_num["dep-CRYCosmicMix"] = [dep_cosmic_run,dep_cosmic_subrun,dep_cosmic_event]
weights = {signal: mean_reco_events[signal]/float(total_sample_events) for signal in norms}
subrun = 0
while True:
	
  events_this_run = 5
  if event_num + events_this_run > total_sample_events:
    events_this_run = total_sample_events - event_num

  datasets = ""
  for signal in norms:
    datasets += "      \"%s\": {\n" % (signal)
    datasets += "        filenames : [\"%s\"]\n" % (filenames[signal])
    datasets += "        weight : %e\n" % (weights[signal])
    if starting_event_num[signal] != [0,0,0]:
      datasets += "        skipToEvent : \"%d:%d:%d\"\n" % (starting_event_num[signal][0],starting_event_num[signal][1],starting_event_num[signal][2])
    datasets += "      }\n"

  d = {}
  d["datasets"] = datasets
  d["outnameMC"] = dirname + "SamplingInput_sr%d_MC.art" % (subrun)
  d["outnameData"] = dirname + "SamplingInput_sr%d_Data.art" % (subrun)

  fout = open(dirname + "SamplingInput_sr%d.fcl" % (subrun),"w")
  fout.write(t.substitute(d))
  fout.close()
  import pdb;pdb.set_trace()
  cmd = ["mu2e","-c",dirname + "SamplingInput_sr%d.fcl" % (subrun),"--nevts","%d" % (events_this_run)]
  p = subprocess.Popen(cmd,stdout=subprocess.PIPE)
  ready = False
  for line in p.stdout:
#    print line
    if "Dataset    Counts    fraction | fraction     weight   Next event" in line:
      ready = True
    if ready:
      if len(line.split()) > 1:
        if line.split()[0].rstrip() in signals:
	  new_run = int(line.rstrip().split()[-5])
	  new_subrun = int(line.rstrip().split()[-3])
	  new_event = int(line.rstrip().split()[-1])
          starting_event_num[line.split()[0]] = [new_run,new_subrun,new_event]
#  print starting_event_num
  p.wait()
  event_num += events_this_run
  print "Job done, return code:",p.returncode,"processed",event_num,"events out of",total_sample_events
  #import pdb;pdb.set_trace()
  if event_num >= total_sample_events:
    break
  subrun+=1

