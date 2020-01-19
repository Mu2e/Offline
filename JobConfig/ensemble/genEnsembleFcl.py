#!/usr/bin/env python
from __future__ import print_function

from string import Template
import sys
import random
import os
from normalizations import *
from argparse import ArgumentParser

def generate(verbose=True,dem_emin=93,dep_emin=83,rpc_tmin=400):
    # when we run from SConscript, the cwd is the python subdir
    # but all file name are relative to Offline, so go there
    cwd = os.getcwd()
    words = cwd.split("/")
    if words[-1] == "ensemble" :
        os.chdir("../..")

    # lists of files to send to scons for dependencies
    sourceFiles = [
	"JobConfig/ensemble/epilog.fcl","JobConfig/ensemble/prolog.fcl",
	"JobConfig/ensemble/epilog_reco.fcl","JobConfig/ensemble/prolog_reco.fcl","JobConfig/ensemble/reco-mcdigis-trig.fcl"]

    targetFiles = []

    projectDir = "gen/fcl/JobConfig/ensemble"
    if not os.path.exists(projectDir) :
      os.makedirs(projectDir)
    
    for tname in ["DIOLeadingLog-cut-mix"]:
      templateFileName = "JobConfig/ensemble/" + tname + ".fcl"
      sourceFiles.append(templateFileName)
      fin = open(templateFileName) 
      t = Template(fin.read())
      d = {"minE": dem_emin, "particleTypes": [11, 13], "minMom": dem_emin}
      fclFileName = projectDir + "/" + tname + ".fcl"
      if verbose:
        print("Creating " + fclFileName)
      targetFiles.append(fclFileName)
      fout = open(fclFileName,"w")
      fout.write(t.substitute(d))
      fout.close()

      templateFileName = "JobConfig/ensemble/reco-mcdigis-trig.fcl"
      fin = open(templateFileName)
      t = Template(fin.read())
      d = {"name": tname}
      fclFileName = projectDir + "/reco-" + tname + ".fcl"
      if verbose:
        print("Creating " + fclFileName)
      targetFiles.append(fclFileName)
      fout = open(fclFileName,"w")
      fout.write(t.substitute(d))
      fout.close()
    
    for tname in ["RPCexternal-cut-mix","RPCinternal-cut-mix"]:
      templateFileName = "JobConfig/ensemble/" + tname + ".fcl"
      sourceFiles.append(templateFileName)
      fin = open(templateFileName) 
      t = Template(fin.read())
      d = {"minE": dep_emin+1, "particleTypes": [11, 13,-11,-13], "minMom": dep_emin, "pionTMin": rpc_tmin}
      fclFileName = projectDir + "/" + tname + ".fcl"
      if verbose:
        print("Creating " + fclFileName)
      targetFiles.append(fclFileName)
      fout = open(fclFileName,"w")
      fout.write(t.substitute(d))
      fout.close()

      templateFileName = "JobConfig/ensemble/reco-mcdigis-trig.fcl"
      fin = open(templateFileName)
      t = Template(fin.read())
      d = {"name": tname}
      fclFileName = projectDir + "/reco-" + tname + ".fcl"
      if verbose:
        print("Creating " + fclFileName)
      targetFiles.append(fclFileName)
      fout = open(fclFileName,"w")
      fout.write(t.substitute(d))
      fout.close()
 
    
    for tname in ["RMCexternal-cut-mix","RMCinternal-cut-mix"]:
      for ikmax in range(8):
        templateFileName = "JobConfig/ensemble/" + tname + ".fcl"
        if ikmax == 0:
          sourceFiles.append(templateFileName)
        fin = open(templateFileName) 
        t = Template(fin.read())
        temp_tname = tname.split("-")[0] + "-kMax%d-" % (ikmax) + tname[len(tname.split("-")[0])+1:] 
    
        d = {"minE": dep_emin+1, "particleTypes": [11, 13,-11,-13], "minMom": dep_emin, "kMaxNum": ikmax, "pionTMin": rpc_tmin}
        fclFileName = projectDir + "/" + temp_tname + ".fcl"
        if verbose:
          print("Creating " + fclFileName)
        targetFiles.append(fclFileName)
        fout = open(fclFileName,"w")
        fout.write(t.substitute(d))
        fout.close()

    for tname in ["RMCexternal-cut-mix","RMCinternal-cut-mix"]:
      for ikmax in range(8):
        temp_tname = tname.split("-")[0] + "-kMax%d-" % (ikmax) + tname[len(tname.split("-")[0])+1:] 
        templateFileName = "JobConfig/ensemble/reco-mcdigis-trig.fcl"
        fin = open(templateFileName)
        t = Template(fin.read())
        d = {"name": temp_tname}
        fclFileName = projectDir + "/reco-" + temp_tname + ".fcl"
        if verbose:
          print("Creating " + fclFileName)
        targetFiles.append(fclFileName)
        fout = open(fclFileName,"w")
        fout.write(t.substitute(d))
        fout.close()

    for tname in ["CeMLeadingLog-mix", "CePLeadingLog-mix"]:
      templateFileName = "JobConfig/ensemble/reco-mcdigis-trig.fcl"
      fin = open(templateFileName)
      t = Template(fin.read())
      d = {"name": tname}
      fclFileName = projectDir + "/reco-" + tname + ".fcl"
      if verbose:
        print("Creating " + fclFileName)
      targetFiles.append(fclFileName)
      fout = open(fclFileName,"w")
      fout.write(t.substitute(d))
      fout.close()
 
    return sourceFiles, targetFiles


if __name__ == "__main__":

    parser = ArgumentParser()
    parser.add_argument("-q", "--quiet",
                        action="store_false", dest="verbose", default=True,
                        help="don't print status messages to stdout")
    parser.add_argument("-t", "--rpc-tmin", dest="rpc_tmin", default=400,
                        help="Early time cutoff for RPC generator")
    parser.add_argument("-p", "--dep-emin", dest="dep_emin", default=83,
			help="Minimum generated momentum for positrons")
    parser.add_argument("-m", "--dem-emin", dest="dem_emin", default=93,
			help="Minimum generated momentum for electrons")
    
    args = parser.parse_args()
    generate(args.verbose)


    exit(0)
      
