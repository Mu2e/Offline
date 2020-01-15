from __future__ import print_function
from string import Template
import sys
import random
import os
from normalizations import *

random.seed()

if len(sys.argv) < 6:
  print("python JobConfig/ensemble/genEnsemble.py <dirname> <max livetime> <livetime for DIO/RPC gen> <kmax number> <run number>")
  sys.exit()

dirname = sys.argv[1]
max_livetime_rmc = float(sys.argv[2]) # in years
max_livetime_others = float(sys.argv[3]) # in years
kmax_number = int(sys.argv[4])
run_number = int(sys.argv[5])

if os.path.exists(os.path.join(os.getcwd(), dirname)):
  print("Error: this directory exists!")
  sys.exit()

os.system("mkdir " + dirname)

livetime_fraction_min = 0.9
livetime_fraction_max = 1.0
livetime = max_livetime_rmc * random.uniform(livetime_fraction_min,livetime_fraction_max)

if max_livetime_rmc <= 6/365.*20/24.:
  print("One week")
  # for one week
  rue_exp_min = -14
  rue_exp_max = -12.8
  rup_exp_min = -14
  rup_exp_max = -12.8
else:
  print("One month")
  # for one month
  rue_exp_min = -14.6
  rue_exp_max = -13.4
  rup_exp_min = -14.6
  rup_exp_max = -13.4


rue = 10**random.uniform(rue_exp_min,rue_exp_max)
rup = 10**random.uniform(rup_exp_min,rup_exp_max)
if rue > 2e-13:
  print("rue too high, change normalization")
  sys.exit()

dem_emin = 93 # generated momentum minimum
dep_emin = 83

tmin = 400 # pion min time for generator

kmax_min = 89
kmax_max = 91

kmax = random.uniform(kmax_min,kmax_max)
if kmax > 91:
  print("kmax too high, change normalization")
  sys.exit()

fout = open(dirname + "/kmax","w")
fout.write("%f\n" % kmax)
fout.close()
fout = open(dirname + "/livetime","w")
fout.write("%f\n" % (livetime*365*24*60*60))
fout.close()
fout = open(dirname + "/rue","w")
fout.write("%e\n" % rue)
fout.close()
fout = open(dirname + "/rup","w")
fout.write("%e\n" % rup)
fout.close()
fout = open(dirname + "/settings","w")
fout.write("%f\n%f\n%f\n%f\n%d\n%d\n" % (dem_emin,dep_emin,tmin,max_livetime_rmc*0.95*365*24*60*60,run_number,random.randint(0,10000)))
fout.close()

# maximum expected events per year
norms = {
  "DIOLeadingLog-cut-mix": dio_normalization(1,dem_emin),
  "CeMLeadingLog-mix": ce_normalization(1,10**rue_exp_max),
  "CePLeadingLog-mix": ce_normalization(1,10**rup_exp_max),
  "RMCexternal-cut-mix": rmc_normalization(1, dep_emin, kmax_max, False),
  "RMCinternal-cut-mix": rmc_normalization(1, dep_emin, kmax_max, True),
  "RPCexternal-cut-mix": 1.59222825e+08, #FIXME python takes too long
  "RPCinternal-cut-mix": 1.098685e+06,
  }

# these have been optimized for 93 MeV and 83 MeV for dem and dep respectively
per_run = {
  "DIOLeadingLog-cut-mix": 250,
  "CeMLeadingLog-mix": 250,
  "CePLeadingLog-mix": 250, 
  "RMCexternal-cut-mix": 300000,
  "RMCinternal-cut-mix": 20000, # for kMax = 91
  "RPCexternal-cut-mix": 300000, 
  "RPCinternal-cut-mix": 2000,

  "reco-DIOLeadingLog-cut-mix": 25,
  "reco-CeMLeadingLog-mix": 25,
  "reco-CePLeadingLog-mix": 25,
  "reco-RMCexternal-cut-mix": 150,
  "reco-RMCinternal-cut-mix": 25,
  "reco-RPCexternal-cut-mix": 10,
  "reco-RPCinternal-cut-mix": 10,
  }

for tname in ["CeMLeadingLog-mix","CePLeadingLog-mix"]:
  fin = open("JobConfig/ensemble/generate_template.sh")
  temp_tname = tname[:-4] + "Mix"
  t = Template(fin.read())

  njobs = int(norms[tname]*max_livetime_others/per_run[tname])+1
  
  d = {"includeOrEmbed": "--include JobConfig/mixing/" + temp_tname + ".fcl", "dirname": dirname, "name": tname, "njobs": njobs, "perjob": per_run[tname]}
  fout = open(dirname + "/generate_" + tname + ".sh","w")
  fout.write(t.substitute(d))
  fout.close()

for tname in ["DIOLeadingLog-cut-mix","RPCexternal-cut-mix","RPCinternal-cut-mix"]:
  fin = open("JobConfig/ensemble/generate_template.sh")
  t = Template(fin.read())

  njobs = int(norms[tname]*max_livetime_others/per_run[tname])+1
  
  d = {"includeOrEmbed": "--embed gen/fcl/JobConfig/ensemble/" + tname + ".fcl", "dirname": dirname, "name": tname, "njobs": njobs, "perjob": per_run[tname]}
  fout = open(dirname + "/generate_" + tname + ".sh","w")
  fout.write(t.substitute(d))
  fout.close()

for tname in ["RMCexternal-cut-mix","RMCinternal-cut-mix"]:
  temp_tname = tname.split("-")[0] + "-kMax%d-" % (kmax_number) + tname[len(tname.split("-")[0])+1:] 
  fin = open("gen/fcl/JobConfig/ensemble/" + temp_tname + ".fcl")
  fout = open(dirname + "/" + tname + ".fcl","w")
  for line in fin:
    fout.write(line)
  fout.write("physics.producers.generate.physics.kMaxUser : %f\n" % kmax)
  fout.close()
  fin.close()
  fin = open("JobConfig/ensemble/generate_template.sh")
  t = Template(fin.read())

  njobs = int(norms[tname]*max_livetime_rmc/per_run[tname])+1
    
  
  d = {"includeOrEmbed": "--embed " + dirname + "/" + tname + ".fcl", "dirname": dirname, "name": temp_tname, "njobs": njobs, "perjob": per_run[tname], "inputs": temp_tname+".txt"}
  fout = open(dirname + "/generate_" + tname + ".sh","w")
  fout.write(t.substitute(d))
  fout.close()


for tname in ["reco-DIOLeadingLog-cut-mix","reco-CeMLeadingLog-mix","reco-CePLeadingLog-mix","reco-RPCexternal-cut-mix","reco-RPCinternal-cut-mix"]:
  fin = open("JobConfig/ensemble/generate_template_reco.sh")
  t = Template(fin.read())
  d = {"includeOrEmbed": "--embed " + "gen/fcl/JobConfig/ensemble/" + tname + ".fcl", "dirname": dirname, "name": tname, "mergeFactor": per_run[tname], "inputs": tname+".txt"}
  fout = open(dirname + "/generate_" + tname + ".sh","w")
  fout.write(t.substitute(d))
  fout.close()

for tname in ["reco-RMCexternal-cut-mix","reco-RMCinternal-cut-mix"]:
  temp_tname = "reco-" + tname.split("-")[1] + "-kMax%d-cut-mix" % (kmax_number)
  fin = open("JobConfig/ensemble/generate_template_reco.sh")
  t = Template(fin.read())
  d = {"includeOrEmbed": "--embed " + "gen/fcl/JobConfig/ensemble/" + temp_tname + ".fcl", "dirname": dirname, "name": temp_tname, "mergeFactor": per_run[tname], "inputs": temp_tname+".txt"}
  fout = open(dirname + "/generate_" + tname + ".sh","w")
  fout.write(t.substitute(d))
  fout.close()


