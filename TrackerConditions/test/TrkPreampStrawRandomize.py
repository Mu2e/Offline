# generates TrkPreampStraw parameters (time delays, thresholds, and adc gains)

import sys
import numpy as np

if len(sys.argv) < 7:
  print("Usage: python TrkPreampStrawRandomize.py <output filename> <time delay sigma (ns)> <mean threshold (mv)> <threshold sigma (mv)> <mean adc gain> <adc gain sigma>")
  sys.exit()

fn = sys.argv[1]

delay_sigma = float(sys.argv[2])
mean_threshold = float(sys.argv[3])
threshold_sigma = float(sys.argv[4])
mean_adc = float(sys.argv[5])
adc_sigma = float(sys.argv[6])

fout = open(fn,"w")
fout.write("# Randomized calibration from TrkPreampStrawRandomize.py, parameters: delay sigma = %f, mean threshold = %f, threshold sigma = %f, mean adc gain = %f, adc gain sigma = %f\n" %
    (delay_sigma, mean_threshold, threshold_sigma, mean_adc, adc_sigma))
fout.write("# python ")
for i in range(len(sys.argv)):
  fout.write(sys.argv[i] + " ")
fout.write("\n#\n#\n#\n")

fout.write("TABLE TrkPreampStraw\n")
fout.write("# index,delay_hv,delay_cal,threshold_hv,threshold_cal,gain\n")
for index in range(36*6*96):
  delay_hv = np.random.normal(0,delay_sigma)
  delay_cal = np.random.normal(0,delay_sigma)
  threshold_hv = np.random.normal(mean_threshold,threshold_sigma)
  threshold_cal = np.random.normal(mean_threshold,threshold_sigma)
  gain = np.random.normal(mean_adc,adc_sigma)
  fout.write("%d, %f, %f, %f, %f, %f\n" % (index, delay_hv, delay_cal, threshold_hv, threshold_cal, gain))

fout.close()
