import random
import sys

if len(sys.argv) != 5:
  print("usage \"randomPhotonYieldDeviation nChannels sigma limitlow limithigh\"")
  print("nChannels=22016 when using the full CRV\"")
  print("nChannels=2048 when using the extracted position\"")
  sys.exit()

nChannels=int(sys.argv[1])
sigma=float(sys.argv[2])
limitlow=float(sys.argv[3])
limithigh=float(sys.argv[4])

print("TABLE CRVPhoton")
print("#{} channels, sigma: {}, limitlow: {}, limithigh: {}".format(nChannels,sigma,limitlow,limithigh))
print("#channel,photonYieldDeviation")

for channel in range(nChannels):
  while True:
    r=random.gauss(0, sigma)
    if r>limitlow and r<limithigh:
      print("{},{:.3f}".format(channel,r))
      break
