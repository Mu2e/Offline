import ROOT
import math
import os

mean_PBI = 3.9e7 # in JobConfig/mixing/prolog.fcl  protonBunchIntensity.extendedMean
dutyfactor = 0.25 # 43.1ms+5ms on spill x8, then 1020ms off spill
ub_per_year = 365*24*60*60./1695e-9*dutyfactor
POT_per_year = ub_per_year*mean_PBI
stopped_mu_per_POT = 0.00187

def cry_normalization(livetime):
  cry_expected_rate = 253440 #Hz
  cry_tmin = 450e-9
  cry_tmax = 1705e-9
  cry_expected_per_ub = cry_expected_rate*(cry_tmax-cry_tmin)
  cry_expected_per_year = ub_per_year * cry_expected_per_ub
  return cry_expected_per_year * livetime

def dayabay_normalization(livetime):
  db_expected_rate = 5494.24 #Hz
  db_tmin = 400e-9
  db_tmax = 1705e-9
  db_expected_per_ub = db_expected_rate*(db_tmax-db_tmin)
  db_expected_per_year = ub_per_year * db_expected_per_ub
  return db_expected_per_year * livetime



def ce_normalization(livetime, rue):
  captures_per_stopped_muon = 0.609
  return POT_per_year * stopped_mu_per_POT * captures_per_stopped_muon * livetime * rue

def dio_normalization(livetime, emin):
  # calculate fraction of spectrum being generated
  spec = open(os.path.join(os.environ["MU2E_BASE_RELEASE"],"ConditionsService/data/czarnecki_szafron_Al_2016.tbl"))
  energy = []
  val = []
  for line in spec:
    energy.append(float(line.split()[0]))
    val.append(float(line.split()[1]))

  total_norm = 0
  cut_norm = 0
  for i in range(len(val)):
    total_norm += val[i]
    if energy[i] >= emin:
      cut_norm += val[i]

  DIO_per_stopped_muon = 0.391 # 1 - captures_per_stopped_muon

  physics_events = POT_per_year * stopped_mu_per_POT * DIO_per_stopped_muon * livetime
  return physics_events * cut_norm/total_norm

def rpc_normalization(livetime, emin, tmin, internal):
  # calculate fraction of spectrum being generated
  spec = open("JobConfig/ensemble/RPCspectrum.dat") # Bistirlich spectrum from 0.05 to 139.95 in steps of 0.1
  energy = []
  val = []
  for line in spec:
    energy.append(float(line.split()[0]))
    val.append(float(line.split()[1]))
  bin_width = energy[1]-energy[0];

  total_norm = 0
  cut_norm = 0
  for i in range(len(val)):
    total_norm += val[i]
    if (energy[i]-bin_width/2. >= emin):
      cut_norm += val[i]

  geometric_stopped_pion_per_POT = 0.002484 # stops assuming infinite pion lifetime (for cd3 sim sample. in docdb-7468)
  RPC_per_stopped_pion = 0.0215; # from reference, uploaded on docdb-469
  internalRPC_per_RPC = 0.00690; # from reference, uploaded on docdb-717


  # calculate survival probability for tmin including smearing of POT
  pot = open("JobConfig/ensemble/POTspectrum.dat"); # ConditionsService/data/potTimingDistribution_20160511.txt, sampled by GenerateProtonTimes_module.cc 
  time = []
  cdf = []
  for line in pot:
    time.append(float(line.split()[0]))
    cdf.append(float(line.split()[1]))
  for i in range(len(cdf)-2,-1,-1):
    cdf[i] += cdf[i+1]
  for i in range(len(cdf)-1,-1,-1):
    cdf[i] /= cdf[0]

  f = ROOT.TFile("/cvmfs/mu2e.opensciencegrid.org/DataFiles/mergedMuonStops/nts.mu2e.pion-DS-TGTstops.MDC2018a.001002_00000000.root");
  d = f.Get("stoppedPionDumper");
  t = d.Get("stops");
  total = 0;
  for i in range(t.GetEntries()):
    #if i%10000 == 0:
    #  print i,t.GetEntries(),i/float(t.GetEntries())
    t.GetEntry(i)
    index = int(tmin - t.time-time[0])
    if (index < 0):
      total += math.exp(-t.tauNormalized);
    elif (index < len(time)-1):
      total += math.exp(-t.tauNormalized)*cdf[index];
  avg_survival_prob = total/t.GetEntries();
  
  physics_events = POT_per_year * geometric_stopped_pion_per_POT * RPC_per_stopped_pion * livetime * avg_survival_prob
  gen_events = physics_events * cut_norm/total_norm;

  if internal:
    physics_events *= internalRPC_per_RPC;
    gen_events *= internalRPC_per_RPC;

  return gen_events

def rmc_normalization(livetime, emin, kmax,internal):
  energy = []
  val = []
  # closure approximation as implemented in MuonCaptureSpectrum.cc
  for i in range(int((kmax-57.05)/0.1)):
    temp_e = 57.05 + i*0.1
    xFit = temp_e/kmax
    energy.append(temp_e)
    val.append((1 - 2*xFit +2*xFit*xFit)*xFit*(1-xFit)*(1-xFit))
  bin_width = energy[1]-energy[0];

  total_norm = 0
  cut_norm = 0
  for i in range(len(val)):
    total_norm += val[i]
    if (energy[i]-bin_width/2. >= emin):
      cut_norm += val[i]

  captures_per_stopped_muon = 0.609
  RMC_gt_57_per_capture = 1.43e-5 # from literature (to overall captures)
  internalRPC_per_RPC = 0.00690; # just copy RPC value

  physics_events = POT_per_year * stopped_mu_per_POT * captures_per_stopped_muon * RMC_gt_57_per_capture * livetime
  gen_events = physics_events * cut_norm/total_norm

  if internal:
    physics_events *= internalRPC_per_RPC;
    gen_events *= internalRPC_per_RPC;

  return gen_events
