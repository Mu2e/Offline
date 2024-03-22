import random
import sys

if len(sys.argv) != 3:
  print("usage \"randomDeadChannels nChannels probability\"")
  print("nChannels=22016 when using the full CRV\"")
  sys.exit()

nChannels=int(sys.argv[1])
probability=float(sys.argv[2])

print("TABLE CRVBadChan")
print("#status 1: not connected")
print("#status 2: dead")
print("#status 3: no data")
print("#status 4: no pedestal")
print("#status 5: no calibration constant")
print("#status 6: noisy")
print("#channel,status")
for channel in range(nChannels):
  if random.random()<probability:
     print("{},2".format(channel))
