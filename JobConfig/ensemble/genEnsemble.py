from string import Template
import sys
import random
import os
from normalizations import *

dirname = sys.argv[1]
max_livetime_rmc = float(sys.argv[2]) # in years
max_livetime_others = float(sys.argv[3]) # in years
kmax_number = int(sys.argv[4])

if os.path.exists(os.path.join(os.getcwd(), dirname)):
  print "Error: this directory exists!"
  sys.exit()

os.system("mkdir " + dirname)
fout = open(dirname + "/settings","w")

livetime_fraction_min = 0.8
livetime_fraction_max = 1.0
livetime = max_livetime_rmc * random.uniform(livetime_fraction_min,livetime_fraction_max)

# for one week
rue_exp_min = -14
rue_exp_max = -12.7
rup_exp_min = -12
rup_exp_max = -10.7

# # for one month
# rue_exp_min = -15
# rue_exp_max = -13.5
# rup_exp_min = -13
# rup_exp_max = -11.5

rue = 10**random.uniform(rue_exp_min,rue_exp_max)
rup = 10**random.uniform(rup_exp_min,rup_exp_max)
if rue > 2e-13:
  print "rue too high, change normalization"
  sys.exit()

dem_emin = 93 # generated momentum minimum
dep_emin = 83

tmin = 400 # pion min time for generator

kmax_min = 89
kmax_max = 91

kmax = random.uniform(kmax_min,kmax_max)
if kmax > 91:
  print "kmax too high, change normalization"
  sys.exit()

fout.write("%f\n%e\n%e\n%f\n%f\n%f\n%f\n" % (livetime,rue,rup,dem_emin,dep_emin,tmin,kmax))
fout.close()

fout = open(dirname + "/kMax%d.fcl" % kmax_number,"w")
fout.close()

# maximum expected events per year
norms = {
  "DIOLeadingLog-cut-mix": dio_normalization(1,dem_emin),
  "CeMLeadingLog-mix": ce_normalization(1,10**rue_exp_max),
  "CePLeadingLog-mix": ce_normalization(1,10**rup_exp_max),
  "RMCexternal-cut": rmc_normalization(1, dep_emin, kmax_max, False),
  "RMCinternal-cut-mix": rmc_normalization(1, dep_emin, kmax_max, True),
  "RPCexternal-cut-mix": 1.59222825e+08, #FIXME python takes too long
  "RPCinternal-cut-mix": 1.098685+06
  }

# these have been optimized for 93 MeV and 83 MeV for dem and dep respectively
per_run = {
  "DIOLeadingLog-cut-mix": 250,
  "CeMLeadingLog-mix": 250,
  "CePLeadingLog-mix": 250, 
  "RMCexternal-cut": 600000,
  "RMCinternal-cut-mix": 20000, # for kMax = 91
  "RPCexternal-cut-mix": 300000, 
  "RPCinternal-cut-mix": 10000
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
  
  d = {"includeOrEmbed": "--include gen/fcl/JobConfig/ensemble/" + tname + ".fcl", "dirname": dirname, "name": tname, "njobs": njobs, "perjob": per_run[tname]}
  fout = open(dirname + "/generate_" + tname + ".sh","w")
  fout.write(t.substitute(d))
  fout.close()

for tname in ["RMCexternal-cut","RMCinternal-cut-mix"]:
  temp_tname = tname.split("-")[0] + "-kMax%d-" % (kmax_number) + tname[len(tname.split("-")[0])+1:] 
  fin = open("gen/fcl/JobConfig/ensemble/" + temp_tname + ".fcl")
  fout = open(dirname + "/" + tname + ".fcl","w")
  for line in fin:
    fout.write(line)
  fout.write("physics.producers.generate.physics.kMaxUser : %f\n" % kmax)
  fout.close()
  fin.close()
  if tname == "RMCexternal-cut":
    fin = open("JobConfig/ensemble/generate_template_nomix.sh")
  else:
    fin = open("JobConfig/ensemble/generate_template.sh")
  t = Template(fin.read())

  njobs = int(norms[tname]*max_livetime_rmc/per_run[tname])+1
    
  
  d = {"includeOrEmbed": "--embed " + dirname + "/" + tname + ".fcl", "dirname": dirname, "name": temp_tname, "njobs": njobs, "perjob": per_run[tname]}
  fout = open(dirname + "/generate_" + tname + ".sh","w")
  fout.write(t.substitute(d))
  fout.close()

