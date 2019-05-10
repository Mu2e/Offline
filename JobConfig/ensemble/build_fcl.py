from string import Template
import sys
import random
import os
from normalizations import *

dirname = sys.argv[1]
max_livetime = float(sys.argv[2]) # in years

if os.path.exists(os.path.join(os.getcwd(), dirname)):
  print "Error: this directory exists!"
  sys.exit()

os.system("mkdir " + dirname)
fout = open(dirname + "/settings","w")

livetime_fraction_min = 0.8
livetime_fraction_max = 1.0
livetime = max_livetime * random.uniform(livetime_fraction_min,livetime_fraction_max)

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

# maximum expected events per year
norms = {
  "dem-DIOLeadingLogMix": dio_normalization(1,dem_emin),
  "dep-DIOLeadingLogMix": 0,
  "dem-CeMLeadingLogMix": ce_normalization(1,10**rue_exp_max),
  "dep-CeMLeadingLogMix": 0,
  "dem-CePLeadingLogMix": 0, 
  "dep-CePLeadingLogMix": ce_normalization(1,10**rup_exp_max),
  "dem-RMCexternalMix": rmc_normalization(1, dem_emin, kmax_max, False),
  "dep-RMCexternalMix": rmc_normalization(1, dep_emin, kmax_max, False),
  "dem-RMCinternalMix": rmc_normalization(1, dem_emin, kmax_max, True),
  "dep-RMCinternalMix": rmc_normalization(1, dep_emin, kmax_max, True),
  "dem-RPCexternalMix": 1.39244803e+08, 
  "dep-RPCexternalMix": 1.59222825e+08,
  "dem-RPCinternalMix": 9.60742e+05,
  "dep-RPCinternalMix": 1.098685+06
  }

# these have been optimized for 93 MeV and 83 MeV for dem and dep respectively
per_run = {
  "dem-DIOLeadingLogMix": 250,
  "dep-DIOLeadingLogMix": 250,
  "dem-CeMLeadingLogMix": 250,
  "dep-CeMLeadingLogMix": 250,
  "dem-CePLeadingLogMix": 250, 
  "dep-CePLeadingLogMix": 250,
  "dem-RMCexternalMix": 400000,
  "dep-RMCexternalMix": 400000,
  "dem-RMCinternalMix": 20000, # for kMax = 91
  "dep-RMCinternalMix": 20000, # for kMax = 91
  "dem-RPCexternalMix": 300000, 
  "dep-RPCexternalMix": 300000,
  "dem-RPCinternalMix": 10000,
  "dep-RPCinternalMix": 10000
  }

for tname in ["CeMLeadingLogMix","DIOLeadingLogMix"]:
  fin = open("JobConfig/ensemble/" + tname + ".fcl") 
  t = Template(fin.read())

  d = {"minE": dem_emin, "particleTypes": [11, 13], "minMom": dem_emin, "dirname": dirname, "name": tname, "part": "dem"}

  fout = open(dirname + "/dem-" + tname + ".fcl","w")
  fout.write(t.substitute(d))
  fout.close()

  fin = open("JobConfig/ensemble/generate_template.sh")
  t = Template(fin.read())

  njobs = int(norms["dem-"+tname]*livetime/per_run["dem-"+tname])+1
  
  d = {"dirname": dirname, "name": tname, "part": "dem", "njobs": njobs, "perjob": per_run["dem-"+tname]}
  fout = open(dirname + "/generate_dem-" + tname + ".sh","w")
  fout.write(t.substitute(d))
  fout.close()

for tname in ["CePLeadingLogMix"]:
  fin = open("JobConfig/ensemble/" + tname + ".fcl") 
  t = Template(fin.read())

  d = {"minE": dep_emin, "particleTypes": [-11, -13], "minMom": dep_emin, "dirname": dirname, "name": tname, "part": "dep"}

  fout = open(dirname + "/dep-" + tname + ".fcl","w")
  fout.write(t.substitute(d))
  fout.close()

  fin = open("JobConfig/ensemble/generate_template.sh")
  t = Template(fin.read())

  njobs = int(norms["dep-"+tname]*livetime/per_run["dep-"+tname])+1
  
  d = {"dirname": dirname, "name": tname, "part": "dep", "njobs": njobs, "perjob": per_run["dep-"+tname]}
  fout = open(dirname + "/generate_dep-" + tname + ".sh","w")
  fout.write(t.substitute(d))
  fout.close()



for tname in ["RPCexternalMix","RPCinternalMix","RMCexternalMix","RMCinternalMix"]:
  fin = open("JobConfig/ensemble/" + tname + ".fcl") 
  t = Template(fin.read())
  
  d = {"minE": dem_emin+1, "particleTypes": [11, 13], "minMom": dem_emin, "kMax": kmax, "pionTMin": tmin, "dirname": dirname, "name": tname, "part": "dem"}
  fout = open(dirname + "/dem-" + tname + ".fcl","w")
  fout.write(t.substitute(d))
  fout.close()

  fin = open("JobConfig/ensemble/generate_template.sh")
  t = Template(fin.read())

  njobs = int(norms["dem-"+tname]*livetime/per_run["dem-"+tname])+1
  
  d = {"dirname": dirname, "name": tname, "part": "dem", "njobs": njobs, "perjob": per_run["dem-"+tname]}
  fout = open(dirname + "/generate_dem-" + tname + ".sh","w")
  fout.write(t.substitute(d))
  fout.close()


  fin = open("JobConfig/ensemble/" + tname + ".fcl") 
  t = Template(fin.read())

  d = {"minE": dep_emin+1, "particleTypes": [-11, -13], "minMom": dep_emin, "kMax": kmax, "pionTMin": tmin, "dirname": dirname, "name": tname, "part": "dep"}
  fout = open(dirname + "/dep-" + tname + ".fcl","w")
  fout.write(t.substitute(d))
  fout.close()

  fin = open("JobConfig/ensemble/generate_template.sh")
  t = Template(fin.read())

  njobs = int(norms["dep-"+tname]*livetime/per_run["dep-"+tname])+1
  
  d = {"dirname": dirname, "name": tname, "part": "dep", "njobs": njobs, "perjob": per_run["dep-"+tname]}
  fout = open(dirname + "/generate_dep-" + tname + ".sh","w")
  fout.write(t.substitute(d))
  fout.close()
