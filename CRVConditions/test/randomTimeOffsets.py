import random
import sys
import csv

if len(sys.argv) != 6:
  print("usage \"randomTimeOffsets nChannels path_to_channelMapFile sigma limitlow limithigh\"")
  print("nChannels=22016 when using the full CRV\"")
  print("nChannels=2048 when using the extracted position\"")
  sys.exit()

nChannels=int(sys.argv[1])
sigma=float(sys.argv[3])
limitlow=float(sys.argv[4])
limithigh=float(sys.argv[5])

with open(sys.argv[2], "r") as channelMapFile:
    fileReader = csv.reader(channelMapFile, delimiter="\t")

    # Skip the first row, which is the header
    next(fileReader)

    timeOffsetsFEB = {}
    timeOffsetsChannel = {}
    for row in fileReader:
        (channelStr, ROC, FEBStr, FEBchannel) = row
        channel=int(channelStr)
        FEB=int(FEBStr)
        if not FEB in timeOffsetsFEB:
           while True:
             r = random.gauss(0, sigma)
             if r>limitlow and r<limithigh:
               timeOffsetsFEB[FEB] = r
               break
        timeOffsetsChannel[channel] = timeOffsetsFEB[FEB]

    print("TABLE CRVTime")
    print("#{} channels, sigma: {}, limitlow: {}, limithigh: {}".format(nChannels,sigma,limitlow,limithigh))
    print("#channel,timeOffset")
    for channel in range(nChannels):
        print("{},{:.3f}".format(channel,timeOffsetsChannel.get(channel,0.0)))

    print("TABLE CRVSiPM")
    print("#channel,pedestal,calibPulseHeight,calibPulseArea")
    for channel in range(nChannels):
        print("{},{:.1f},{:.1f},{:.1f}".format(channel,100.0,11.4,394.6))
