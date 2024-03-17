import sys
import os
import ROOT
import numpy as np
import matplotlib.pyplot as plt

rootFile = ROOT.TFile.Open("FPGAratesNew.root","READ")
rootFile.cd("CrvRates/FPGAhitRates")
nFPGAs=ROOT.gDirectory.GetListOfKeys().GetEntries()
nEvents=ROOT.gDirectory.FindObjectAny(ROOT.gDirectory.GetListOfKeys().At(0).GetName()).GetEntries()
print(f"{nFPGAs} {nEvents}")

#create csv file of all FPGAs for hit thresholds up to 60
for iFPGA in range(nFPGAs):
  FPGAname=ROOT.gDirectory.GetListOfKeys().At(iFPGA).GetName()
  print(FPGAname, end='')
  for threshold in range(1,61):
    eventsAboveThreshold=ROOT.gDirectory.FindObjectAny(FPGAname).Integral(threshold,100)
    eventsLost=eventsAboveThreshold/nEvents
    print(f",{eventsLost}", end='')
  print("") #for new line

#create plot for fraction of events above threshold for selected FPGAs
plt.figure(figsize=(12,6))
plt.title("Fraction of on-spill events above FPGA hit threshold (at individual FPGAs)")
plt.xlabel("Maximum number of hits per FPGA and microbunch")
plt.ylabel("Fraction of on-spill events above hit threshold")
plt.xlim(0,60)
plt.ylim(0,1)

selectedFPGAindices=[(809,"-g","Worst CRV-T FPGA"),(751,"-b","Worst CRV-TS FPGA"),(1105,"--c","Worst CRV-TS-Ext FPGA"),(1232,"--k","Worst CRV-U FPGA"),(15,"-m","Worst CRV-R FPGA"),(455,"-y","Worst CRV-L FPGA")]
thresholds = [i for i in range(1,61)]
for iFPGA in selectedFPGAindices:
  FPGAname=ROOT.gDirectory.GetListOfKeys().At(iFPGA[0]).GetName()
  eventsLost = []
  for threshold in range(1,61):
    eventsAboveThreshold=ROOT.gDirectory.FindObjectAny(FPGAname).Integral(threshold,100)
    eventsLost.append(eventsAboveThreshold/nEvents)
  plt.plot(thresholds, eventsLost, iFPGA[1], label=iFPGA[2])

plt.legend(loc="upper right", frameon=False)
plt.text(20,0.7,"uses only 4 channels per FPGA for CRV-U and CRV-TS-Ext")
#plt.show()

#plot FPGA hit rates of selected FPGAs
canvas = ROOT.TCanvas("FPGA hit rates","FPGA hit rates",1200,600)
canvas.Divide(3,2)
for i in range(len(selectedFPGAindices)):
  canvas.cd(i+1)
  ROOT.gPad.SetLogy()
  iFPGA=selectedFPGAindices[i]
  FPGAname=ROOT.gDirectory.GetListOfKeys().At(iFPGA[0]).GetName()
  hist=ROOT.gDirectory.FindObjectAny(FPGAname)
  print(hist.GetTitle())
  hist.SetTitle(iFPGA[2])
  hist.GetXaxis().SetRangeUser(0,60)
  hist.GetXaxis().SetTitle("number of hits per FPGA and microbunch")
  hist.GetYaxis().SetTitle("number of on-spill events")
  hist.Draw()
  ROOT.gPad.Modified()
  ROOT.gPad.Update()
  FPGA4channels=["CRV-TS-Ext","CRV-U"]
  if any(s in iFPGA[2] for s in FPGA4channels):
     t=ROOT.TText(0.3,0.7,"uses only 4 channels of the FPGA")
     t.SetNDC()
     t.SetTextFont(42)
     t.DrawClone()
  else:
     t=ROOT.TText(0.3,0.7,"uses all 16 channels of the FPGA")
     t.SetNDC()
     t.SetTextFont(42)
     t.DrawClone()

#plot FPGA hit multiplicities of selected FPGAs
rootFile.cd("CrvRates/FPGAhitMultiplicities")
canvas2 = ROOT.TCanvas("FPGA hit multiplicities","FPGA hit multiplicities",1200,600)
canvas2.Divide(3,2)
for i in range(len(selectedFPGAindices)):
  canvas2.cd(i+1)
  ROOT.gPad.SetLogy()
  iFPGA=selectedFPGAindices[i]
  FPGAname=ROOT.gDirectory.GetListOfKeys().At(iFPGA[0]).GetName()
  hist=ROOT.gDirectory.FindObjectAny(FPGAname)
  print(hist.GetTitle())
  hist.SetTitle(iFPGA[2])
  hist.GetXaxis().SetRangeUser(0,20)
  hist.GetXaxis().SetTitle("number of subsequent \"hits\" within a pulse")
  hist.Draw()
  ROOT.gPad.Modified()
  ROOT.gPad.Update()

plt.show()
