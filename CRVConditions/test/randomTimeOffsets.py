import random
import sys
import csv

if len(sys.argv) != 3:
  print("usage \"randomTimeOffsets nChannels path_to_channelMapFile\"")
  print("nChannels=22016 when using the full CRV\"")
  print("nChannels=2048 when using the extracted position\"")
  sys.exit()

nChannels=int(sys.argv[1])

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
           timeOffsetsFEB[FEB] = random.gauss(0, 1)
        timeOffsetsChannel[channel] = timeOffsetsFEB[FEB]

    print("TABLE CRVTime")
    print("#channel,timeOffset")
    for channel in range(nChannels):
        print("{},{:.3f}".format(channel,timeOffsetsChannel.get(channel,0.0)))

    print("TABLE CRVSiPM")
    print("#channel,pedestal,calibPulseHeight,calibPulseArea")
    for channel in range(nChannels):
        print("{},{:.1f},{:.1f},{:.1f}".format(channel,100.0,11.4,394.6))
