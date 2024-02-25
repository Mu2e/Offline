import random
import sys
import csv

if len(sys.argv) != 2:
  print("usage \"randomTimeOffsets path_to_channelMapFile\"")
  sys.exit()

with open(sys.argv[1], "r") as channelMapFile:
    fileReader = csv.reader(channelMapFile, delimiter="\t")

    # Skip the first row, which is the header
    next(fileReader)

    print("TABLE CRVTime")
    print("#channel,timeOffset")

    channels = []
    timeOffsets = {}
    for row in fileReader:
        (channel, ROC, FEB, FEBchannel) = row
        channels.append(channel)
        if not FEB in timeOffsets:
           timeOffsets[FEB] = random.gauss(0, 1)
        print("{},{:.3f}".format(channel,timeOffsets[FEB]))

    print("TABLE CRVSiPM")
    print("#channel,pedestal,calibPulseHeight,calibPulseArea")
    for channel in channels:
        print("{},{:.1f},{:.1f},{:.1f}".format(channel,100.0,11.4,394.6))
