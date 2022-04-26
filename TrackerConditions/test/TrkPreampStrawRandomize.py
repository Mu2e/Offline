# generates TrkPreampStraw parameters (time delays, thresholds, and adc gains)

import sys
import numpy as np

NPlanes = 36
NPanels = 6
NStraws = 96
NUStraws = NPlanes*NPanels*NStraws

if len(sys.argv) < 7:
  print("Usage: python TrkPreampStrawRandomize.py <output filename> <seed> <absolute time offset sigma (ns)> <deltat time offset sigma (ns)> <cal mean threshold (mv)> <mean hv-cal threshold (mv)> <cal threshold sigma (mv)> <hv threshold sigma (mv)> <mean adc gain> <adc gain sigma>")
  sys.exit()

fn = sys.argv[1]

seed = int(sys.argv[2])
abs_delay_sigma = float(sys.argv[3])
dt_delay_sigma = float(sys.argv[4])
mean_threshold = float(sys.argv[5])
threshold_offset = float(sys.argv[6])
cal_threshold_sigma = float(sys.argv[7])
hv_threshold_sigma = float(sys.argv[8])
mean_adc = float(sys.argv[9])
adc_sigma = float(sys.argv[10])

fout = open(fn,"w")
fout.write("# Randomized calibration from TrkPreampStrawRandomize.py, parameters:\n")
fout.write("# absolute time offset sigma = %f ns\n" % (abs_delay_sigma))
fout.write("# deltat time offset sigma = %f ns\n" % (dt_delay_sigma))
fout.write("# cal mean threshold = %f mv\n" % (mean_threshold))
fout.write("# mean hv-cal threshold = %f mv\n" % (threshold_offset))
fout.write("# cal threshold sigma = %f mv\n" % (cal_threshold_sigma))
fout.write("# hv threshold sigma = %f mv\n" % (hv_threshold_sigma))
fout.write("# mean adc gain = %f\n" % (mean_adc))
fout.write("# adc gain sigma = %f %%\n" % (adc_sigma))
fout.write("# python ")
for i in range(len(sys.argv)):
  fout.write(sys.argv[i] + " ")
fout.write("\n#\n#\n#\n")

rng = np.random.default_rng(seed)

fout.write("TABLE TrkPreampStraw\n")
fout.write("# index,delay_hv,delay_cal,threshold_hv,threshold_cal,gain\n")
abs_delay = rng.normal(0,abs_delay_sigma,NUStraws)
dt_delay = rng.normal(0,dt_delay_sigma,NUStraws)
threshold_hv = rng.normal(mean_threshold+threshold_offset,hv_threshold_sigma,NUStraws)
threshold_cal = rng.normal(mean_threshold,cal_threshold_sigma,NUStraws)
gain = rng.normal(mean_adc,adc_sigma/100.*mean_adc,NUStraws)
for i in range(NUStraws):
  delay_hv = abs_delay[i] + dt_delay[i]/2.
  delay_cal = abs_delay[i] - dt_delay[i]/2.
  fout.write("%d, %f, %f, %f, %f, %f\n" % (i, delay_hv, delay_cal, threshold_hv[i], threshold_cal[i], gain[i]))

fout.close()
