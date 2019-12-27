from __future__ import print_function
from string import Template
import sys
import random
import os
import glob
import ROOT
from normalizations import *
import subprocess

max_events_per_subrun = 1000000

dirname = sys.argv[1]
outpath = sys.argv[2]
livetime = float(open(os.path.join(dirname,"livetime")).readline())/365./24./60./60.
rue = float(open(os.path.join(dirname,"rue")).readline())
rup = float(open(os.path.join(dirname,"rup")).readline())
kmax = float(open(os.path.join(dirname,"kmax")).readline())

fin = open(os.path.join(dirname,"settings"))
lines = fin.readlines()

dem_emin = float(lines[0])
dep_emin = float(lines[1])
tmin = float(lines[2])
max_livetime = float(lines[3])/365./24./60./60.
run = int(lines[4])
samplingseed = int(lines[5])

#use_dayabay = 0
#if len(sys.argv) > 2:
#  use_dayabay = int(sys.argv[2])

ROOT.gRandom.SetSeed(0)

norms = {
  "DIOLeadingLog": dio_normalization(livetime,dem_emin),
  "CeMLeadingLog": ce_normalization(livetime,rue),
  "CePLeadingLog": ce_normalization(livetime,rup),
  "RMCexternal": rmc_normalization(livetime, dep_emin+1, kmax, False), # this is expected number of gammas
  "RMCinternal": rmc_normalization(livetime, dep_emin+1, kmax, True),
  "RPCexternal": 1.59222825e+08*livetime, #FIXME python takes too long
  "RPCinternal": 1.098685e+06*livetime,
  "DYBCosmic": dayabay_normalization(max_livetime),
  "CRYCosmic": cry_normalization(max_livetime)
  }


starting_event_num = {}
currently_used_events = {}
max_possible_events = {}
mean_reco_events = {}
filenames = {}

ffns = open(os.path.join(dirname,"filenames"))
for line in ffns: 
  signal = line.split()[0]
  print(signal)
  filenames[signal] = line.split()[1]
  starting_fraction = float(line.split()[2])
  max_fraction = float(line.split()[3])
  fin = ROOT.TFile(filenames[signal])
  te = fin.Get("Events")
  
  # determine total number of events surviving all cuts
  reco_events = te.GetEntries()

  # get event id for first event we will use
  start_entry = int(reco_events * starting_fraction)
  tfi = fin.Get("FileIndex")
  current_entry = 0
  check_entry = 0
  for i in range(tfi.GetEntries()):
    tfi.GetEntry(i)
    if tfi.Element.getEntryType() != 2:
      continue
    if current_entry >= start_entry:
      start_run_id = tfi.Element.eventID_.run()
      start_subrun_id = tfi.Element.eventID_.subRun()
      start_event_id = tfi.Element.eventID_.event()
      check_entry = tfi.Element.entry_
      if start_subrun_id == 0:
        # fails if we try to use subrun of 0 so get next event
        continue
      break
    current_entry += 1

  te.GetEntry(check_entry)
  if (not start_run_id == te.EventAuxiliary.run() or not start_subrun_id == te.EventAuxiliary.subRun() or not start_event_id == te.EventAuxiliary.event()):
    print("PROBLEM GETTING STARTING EVENT")
    import pdb;pdb.set_trace()
  starting_event_num[signal] = [start_run_id,start_subrun_id,start_event_id]
  # determine how many events we can use before throwing error
  max_possible_events[signal] = int(reco_events * max_fraction) - int(reco_events * starting_fraction)
  currently_used_events[signal] = 0

  # determine total number of events generated
  t = fin.Get("SubRuns")
  # find the right branch
  bl = t.GetListOfBranches()
  bn = ""
  for i in range(bl.GetEntries()):
    if bl[i].GetName().split("__")[0] == "mu2e::GenEventCount_genCounter":
      bn = bl[i].GetName()
  gen_events = 0
  for i in range(t.GetEntries()):
    t.GetEntry(i)
    gen_events += getattr(t,bn).count()
  
 
  mean_gen_events = norms[signal]
  mean_reco_events[signal] = mean_gen_events*reco_events/float(gen_events)
  # DEBUG ONLY
  # print signal,"GEN_EVENTS:",gen_events,"RECO_EVENTS:",reco_events,"EXPECTED EVENTS:",mean_reco_events[signal],"STARTING AT:",starting_event_num[signal]

total_sample_events = ROOT.gRandom.Poisson(sum(mean_reco_events.values()))
# DEBUG ONLY
# print "TOTAL EXPECTED EVENTS:",sum(mean_reco_events.values()),"GENERATING:",total_sample_events

# calculate the normalized weights for each signal
weights = {signal: mean_reco_events[signal]/float(total_sample_events) for signal in mean_reco_events}

# generate subrun by subrun

fin = open("JobConfig/ensemble/SamplingInput.fcl")
t = Template(fin.read())

subrun = 0

num_events_already_sampled = 0

while True:
	
  events_this_run = max_events_per_subrun
  if num_events_already_sampled + events_this_run > total_sample_events:
    events_this_run = total_sample_events - num_events_already_sampled

  datasets = ""
  for signal in weights:
    datasets += "      %s: {\n" % (signal)
    datasets += "        fileNames : [\"%s\"]\n" % (filenames[signal])
    datasets += "        weight : %e\n" % (weights[signal])
    if starting_event_num[signal] != [0,0,0]:
      datasets += "        skipToEvent : \"%d:%d:%d\"\n" % (starting_event_num[signal][0],starting_event_num[signal][1],starting_event_num[signal][2])
    datasets += "      }\n"

  d = {}
  d["datasets"] = datasets
  d["outnameMC"] = os.path.join(outpath,"mcs.mu2e.ensemble-MC.MDC2018i.%06d_%08d.art" % (run,subrun))
  d["outnameData"] = os.path.join(outpath,"mcs.mu2e.ensemble-Data.MDC2018i.%06d_%08d.art" % (run,subrun))
  d["run"] = run
  d["subRun"] = subrun
  d["samplingSeed"] = samplingseed + subrun
  # put all the exact parameter values in the fcl file
  d["comments"] = "#livetime: %f\n#rue: %e\n#rup: %e\n#kmax: %f\n#dem_emin: %f\n#dep_emin: %f\n#tmin: %f\n#max_livetime: %f\n#run: %d\n" % (livetime*365*24*60*60,rue,rup,kmax,dem_emin,dep_emin,tmin,max_livetime*365*24*60*60,run)

  fout = open(os.path.join(dirname,"SamplingInput_sr%d.fcl" % (subrun)),"w")
  fout.write(t.substitute(d))
  fout.close()

  flog = open(os.path.join(dirname,"SamplingInput_sr%d.log" % (subrun)),"w")
  
  cmd = ["mu2e","-c",os.path.join(dirname,"SamplingInput_sr%d.fcl" % (subrun)),"--nevts","%d" % (events_this_run)]
  p = subprocess.Popen(cmd,stdout=subprocess.PIPE)
  ready = False
  for line in p.stdout:
    flog.write(line)
#    DEBUG only
#    print line
    if "Dataset" in line.split() and "Counts" in line.split() and "fraction" in line.split() and "Next" in line.split():
#    if "Dataset    Counts    fraction | fraction     weight   Next event" in line:
      ready = True
    if ready:
      if len(line.split()) > 1:
        if line.split()[0].rstrip() in starting_event_num:
	  new_run = int(line.rstrip().split()[-5])
	  new_subrun = int(line.rstrip().split()[-3])
	  new_event = int(line.rstrip().split()[-1])
          starting_event_num[line.split()[0]] = [new_run,new_subrun,new_event]
          currently_used_events[line.split()[0]] += int(line.split()[1])
  p.wait()
  num_events_already_sampled += events_this_run
  print("Job done, return code: " + p.returncode + " processed " + num_events_already_sampled + " events out of " + total_sample_events)
  if num_events_already_sampled >= total_sample_events:
    break
  for signal in currently_used_events:
    if currently_used_events[signal] > max_possible_events[signal]:
      print("SIGNAL " + signal + " HAS RUN OUT OF EVENTS! " + currently_used_events[signal] + " " + max_possible_events[signal])
    print("SIGNAL " + signal + " " + currently_used_events[signal] + " " + max_possible_events[signal])

  subrun+=1

