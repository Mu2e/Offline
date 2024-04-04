import random
import sys
import math

if len(sys.argv) != 3:
  print("usage \"randomDeadChannels nChannels probability\"")
  print("nChannels=22016 when using the full CRV\"")
  print("nChannels=2048 when using the extracted position\"")
  sys.exit()

nChannels=int(sys.argv[1])
probability=float(sys.argv[2])

print("TABLE CRVBadChan")
print("#{} channels, probability: {}".format(nChannels,probability))
print("#A channel can have more than one status bit. Add all status bits for the complete status.")
print("#status 1 (bit 0): not connected")
print("#status 2 (bit 1): ignore channel in reconstruction")
print("#status 4 (bit 2): no data")
print("#status 8 (bit 3): no pedestal")
print("#status 16 (bit 4): no calibration constant")
print("#status 32 (bit 5): noisy")
print("#channel,status")
prevDeadCounter=-2
for channel in range(nChannels):
  if random.random()<probability:
     deadCounter=math.floor(channel/4)
     if deadCounter-prevDeadCounter>1:
        prevDeadCounter=deadCounter
        print("{},2".format(channel))
