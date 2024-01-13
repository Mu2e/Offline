import random
import sys

if len(sys.argv) != 2:
  print("usage \"randomScintDeviation nChannels\"")
  print("nChannels=22016 when using the full CRV\"")
  sys.exit()

nChannels=int(sys.argv[1])

print("TABLE CRVScint")
print("#channel,scintYieldDeviation")
for channel in range(nChannels):
  print("{},{:.3f}".format(channel,random.gauss(0, 1)))
