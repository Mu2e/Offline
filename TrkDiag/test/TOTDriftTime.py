import ROOT
import sys
import numpy as np

# usage: python TOTDriftTime.py <name of root file with shdiag ttree> <num energy bins> <energy binning (KeV)> <num tot bins> <TOT binning (ns)>

pcode = int(sys.argv[2])
ebins = int(sys.argv[3])
ewidth = float(sys.argv[4])
tbins = int(sys.argv[5])
twidth = float(sys.argv[6])

f = ROOT.TFile(sys.argv[1])
t = f.Get("SHD/shdiag")

means = np.zeros((tbins,ebins))

hdrift = ROOT.TH1F("hdrift","hdrift",200,-100,100)
t.Project("hdrift","ctime-mcsptime+pbtmc-ptime","mcsptime < 1800 && mcproc == %d" % pcode)
default_drift = hdrift.GetMean()

htn = ROOT.TH1F("htn","htn",tbins,0,tbins*twidth)
htd = ROOT.TH1F("htd","htd",tbins,0,tbins*twidth)
t.Project("htn","totcal*(tcal < thv)+tothv*(tcal >= thv)","(ctime-mcsptime+pbtmc-ptime)*(mcsptime < 1800 && mcproc == %d)" % pcode)
t.Project("htd","totcal*(tcal < thv)+tothv*(tcal >= thv)","(mcsptime < 1800 && mcproc == 56)")
htn.Divide(htd)

hn = ROOT.TH2F("hn","hn",tbins,0,tbins*twidth,ebins,0,ebins*ewidth)
hd = ROOT.TH2F("hd","hd",tbins,0,tbins*twidth,ebins,0,ebins*ewidth)
t.Project("hn","edep:totcal*(tcal < thv)+tothv*(tcal >= thv)","(ctime-mcsptime+pbtmc-ptime)*(mcsptime < 1800 && mcproc == %d)" % pcode)
t.Project("hd","edep:totcal*(tcal < thv)+tothv*(tcal >= thv)","(mcsptime < 1800 && mcproc == %d)" % pcode)
hn.Divide(hd)
for i in range(tbins):
  for j in range(ebins):
    means[i][j] = hn.GetBinContent(i+1,j+1)
    if hd.GetBinContent(i+1,j+1) < 10:
      if htd.GetBinContent(i+1,j+1) < 10:
        means[i][j] = default_drift
      else:
        means[i][j] = htn.GetBinContent(i+1)
    print("%.2f, " % means[i][j],end='')
  print("")

hval = ROOT.TH1F("hval","hval",100,-50,50)
ht2d = ROOT.TH2F("ht2d","ht2d",100,-50,50,tbins,0,tbins*twidth)
he2d = ROOT.TH2F("he2d","he2d",100,-50,50,ebins,0,ebins*ewidth)

for i in range(t.GetEntries()):
  if i > 100000:
    break
  t.GetEntry(i)
  if t.mcsptime > 1800 or t.mcproc != pcode:
    continue
  ebin = max(0,min(ebins-1,int(t.edep/ewidth)))
  if t.tcal < t.thv:
    tbin = int(t.totcal/twidth)
  else:
    tbin = int(t.tothv/twidth)
  tbin = max(0,min(tbins-1,tbin))
  val = means[tbin][ebin]
  hval.Fill((t.ctime-t.mcsptime+t.pbtmc-t.ptime)-val)
  ht2d.Fill((t.ctime-t.mcsptime+t.pbtmc-t.ptime)-val,tbin*twidth)
  he2d.Fill((t.ctime-t.mcsptime+t.pbtmc-t.ptime)-val,ebin*ewidth)
print("PROCESSING DONE")
htp = ht2d.ProfileY("htp",1,-1,"S")
htp.GetXaxis().SetTitle("TOT (ns)")
htp.GetYaxis().SetTitle("Mean drift time residual (ns)")
htp.Draw()
input()
hep = he2d.ProfileY("hep",1,-1,"S")
hep.GetXaxis().SetTitle("edep (KeV)")
hep.GetYaxis().SetTitle("Mean drift time residual (ns)")
hep.Draw()
input()
hval.GetXaxis().SetTitle("Drift time residual (ns)")
hval.Draw()
input()
