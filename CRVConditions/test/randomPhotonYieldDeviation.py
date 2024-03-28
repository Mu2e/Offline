import random
import sys

if len(sys.argv) != 2:
  print("usage \"randomPhotonYieldDeviation nChannels\"")
  print("nChannels=22016 when using the full CRV\"")
  print("nChannels=2048 when using the extracted position\"")
  sys.exit()

nChannels=int(sys.argv[1])

print("TABLE CRVPhoton")
print("#channel,photonYieldDeviation")
for channel in range(nChannels):
  print("{},{:.3f}".format(channel,random.gauss(0, 1)))
