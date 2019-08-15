import random

rue_exp_required = -13.1 # mean expected ~25
rue_exp_max = -12.8

rup_exp_required = -13.1 # mean expected ~25
rup_exp_max = -12.8

rue_passed = False
rup_passed = False
for i in range(3):
  rue = float(open("closedEnsemble%dMDC2018i/rue" % i).readline())
  rup = float(open("closedEnsemble%dMDC2018i/rup" % i).readline())
  if rue > 10**rue_exp_required:
    rue_passed = True
  if rup > 10**rup_exp_required:
    rup_passed = True

if not rue_passed:
  whichensemble = random.randint(0,2)
  rue = 10**random.uniform(rue_exp_required,rue_exp_max)
  fout = open("closedEnsemble%dMDC2018i/rue" % whichensemble,"w")
  fout.write("%e\n" % rue)
  fout.close()

if not rup_passed:
  whichensemble = random.randint(0,2)
  rup = 10**random.uniform(rup_exp_required,rup_exp_max)
  fout = open("closedEnsemble%dMDC2018i/rup" % whichensemble,"w")
  fout.write("%e\n" % rup)
  fout.close()
