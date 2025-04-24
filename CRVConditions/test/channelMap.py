sectors = {
            "R1": 5,
            "R2": 13,
            "R3": 1,
            "R4": 1,
            "R5": 2,
            "R6": 4,
            "L1": 4,
            "L2": 11,
            "L3": 4,
            "T1": 3,
            "T2": 2,
            "T3": 4,
            "T4": 13,
            "T5": 3,
            "E1": 1,
            "E2": 1,
            "U" : 4,
            "D1": 3,
            "D2": 1,
            "D3": 1,
            "D4": 1,
            "C1": 1,
            "C2": 1
          }

#for special / out-of-sequence situations
outOfSequence = {
                  "R3": {"ROC": 4,  "port": 5},     #module above cryo hole
                  "R4": {"ROC": 4,  "port": 1},     #module below cryo hole
                  "R5": {"ROC": 4,  "port": 13},    #modules downstream of cryo hole
                  "R6": {"ROC": 5,  "port": 1},     #endcap R modules
                  "L1": {"ROC": 6,  "port": 1},     #L modules
                  "L3": {"ROC": 9,  "port": 1},     #endcap L modules
                  "T1": {"ROC": 10, "port": 1},     #TS modules
                  "T3": {"ROC": 10, "port": 17},    #T modules
                  "E1": {"ROC": 14, "port": 1},     #E modules
                  "U":  {"ROC": 15, "port": 1},     #U modules
                  "D1": {"ROC": 17, "port": 1},     #D modules
                  "C1": {"ROC": 4,  "port": 9},     #long cryo module
                  "C2": {"ROC": 4,  "port": 12}     #short cryo module
                }

#one-sided readouts (they don't use FEBs at the opposite side, and therefore only need 2 FEBs per module)
oneSidedReadouts = {"R3": 1, "T1": 0, "T2": 0, "E1": 1, "E2": 1, "U": 0, "D2": 0, "D3": 1, "C1": 1, "C2": 1}

#sparsified FEBs
sparsifiedFEBs = {"E1": 4, "E2": 4, "U": 4}  #all with a factor of 4

nCountersPerModule=16*4
nCountersPerModuleC1=20*4
nCountersPerModuleC2=4*4
nSiPMsPerCounter=4

globalChannel=0
ROC=1  #starts counting at 1 - same as CAT6 cable labels
FEB=0
FEBchannel=0

print("Channel\tROC\tFEB\tFEBchannel")
for sectorName,sectorModules in sectors.items():
  #deal with special / out-of-sequence situations
  if sectorName in outOfSequence:
    ROC=outOfSequence[sectorName]["ROC"]
    FEB=outOfSequence[sectorName]["port"]-1  #FEB start counting at 0, ports start counting at 1 (as printed at ROC and CAT6 cable label)

  for sectorModule in range(sectorModules):
    #print("sector {} module {}".format(sectorName,sectorModule))

    #special situation for U modules
    if sectorName=="U" and sectorModule==2:
       ROC=16  #ROCs start counting at 1
       FEB=0   #FEBs start counting at 0, ports start counting at 1 (as printed at ROC and CAT6 cable label)

    #special situation for C module
    nCountersThisModule=nCountersPerModule
    if sectorName=="C1":
       nCountersThisModule=nCountersPerModuleC1
    if sectorName=="C2":
       nCountersThisModule=nCountersPerModuleC2

    for moduleCounter in range(nCountersThisModule):
       for iSiPM in range(nSiPMsPerCounter):

          #skip channels for oneSidedReadouts
          if sectorName in oneSidedReadouts:
            if iSiPM%2!=oneSidedReadouts[sectorName]:
               globalChannel+=1    #global channel still needs to be increased
               continue

          #the SiPMs at the negative counter end get FEB numbers FEB+0 and FEB+1.
          #the SiPMs at the positive counter end get FEB numbers FEB+2 and FEB+3.
          #the above does not apply to one-sided readouts
          actualFEB=FEB
          if iSiPM%2==1 and (sectorName not in oneSidedReadouts):
             actualFEB+=2

          FEBchannel=moduleCounter*2+iSiPM//2
          if sectorName in sparsifiedFEBs:
             FEBchannel*=sparsifiedFEBs[sectorName]
          #FEBs have only 64 channel. go to the 2nd FEB for channel numbers above 64.
          actualFEBchannel=FEBchannel%64
          actualFEB+=FEBchannel//64
          print("{}\t{}\t{}\t{}".format(globalChannel,ROC,actualFEB,actualFEBchannel))
          globalChannel+=1

    if sectorName=="C1":    #cryo modules are special situation
       FEB+=3
    elif sectorName=="C2":
       FEB+=1
    else:  #regular modules
       if sectorName not in oneSidedReadouts:
          if sectorName not in sparsifiedFEBs:
             FEB+=4  #each module requires 4 FEBs. jump ahead by 4 FEBs for the next module.
          else:
             FEB+=4*sparsifiedFEBs[sectorName]
       else:
          if sectorName not in sparsifiedFEBs:
             FEB+=2  #each module requires 2 FEBs for one-sided readout modules. jump ahead by 2 FEBs for the next module.
          else:
             FEB+=2*sparsifiedFEBs[sectorName]

    if FEB>=24:   #each ROC can handle 24 FEB. use the next ROC if all FEBs have been used
       FEB=0
       ROC+=1
