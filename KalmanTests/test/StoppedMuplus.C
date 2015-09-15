//#include "TROOT.h"
//#include "TH1F.h"
//#include "TLegend.h"
//#include "TChain.h"
{
  gROOT->LoadMacro("../../Scripts/FillChain.C");
  TChain* mup(0);
  FillChain(mup,"/data/HD2/brownd/mu2e/StoppedMuplus.3795378","StoppedMuplus.root",100,"RKFDownstreamePlus/trkdiag");
  TCut mupsel("trkqual>0.8&&momerr<0.4");
  TH1F* recmom= new TH1F("recmom","Reconstructed Momentum;Fit Momentum (MeV/c)",50,45,55);
  recmom->SetStats(0);
  TH2F* mpos = new TH2F("mpos","Muon Decay Position of Good Tracks;x (mm); y (mm)",50,-80,80,50,-80,80);
  mpos->SetStats(0);
  mup->Project("recmom","fit.mom",mupsel);
  mup->Project("mpos","mcgen.y:mcgen.x",mupsel);
  TCanvas* can = new TCanvas("momcan","momcan",600,1200);
  can->Divide(1,2);
  can->cd(1);
  mpos->Draw("colorz");
  can->cd(2);
  recmom->Draw();

}
