#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile2D.h"
#include "TLegend.h"
#include "TTree.h"
#include "THStack.h"
#include "TH1F.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TCut.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TLine.h"
#include "TRandom3.h"
#include <iostream>
#include <math.h>
#include <vector>

void StereoPositionTest(TTree* spdiag) {
  gStyle->SetOptStat(0);
  TCut stereo("stereo&&mvaout>0.7");
  TCut de("mcproc<20");
  TCut ce("mcgen==2");
  TCut ldist("dz<40");
  TCut gdist("dz>40");
  TH1F* drstcel = new TH1F("drstcel","Reco - True transverse radius, #Delta z< 40 mm, CE;#Delta #rho (mm)",100,-100,100);
  TH1F* drstdel = new TH1F("drstdel","Reco - True transverse radius, #Delta z< 40 mm, #delta e^{-};#Delta #rho (mm)",100,-100,100);
  TH1F* drshcel = new TH1F("drshcel","Reco - True transverse radius, #Delta z< 40 mm, CE;#Delta #rho (mm)",100,-100,100);
  TH1F* drshdel = new TH1F("drshdel","Reco - True transverse radius, #Delta z< 40 mm, #delta e^{-};#Delta #rho (mm)",100,-100,100);
  TH1F* drstceg = new TH1F("drstceg","Reco - True transverse radius, #Delta z> 40 mm, CE;#Delta #rho (mm)",100,-100,100);
  TH1F* drstdeg = new TH1F("drstdeg","Reco - True transverse radius, #Delta z> 40 mm, #delta e^{-};#Delta #rho (mm)",100,-100,100);
  TH1F* drshceg = new TH1F("drshceg","Reco - True transverse radius, #Delta z> 40 mm, CE;#Delta #rho (mm)",100,-100,100);
  TH1F* drshdeg = new TH1F("drshdeg","Reco - True transverse radius, #Delta z> 40 mm, #delta e^{-};#Delta #rho (mm)",100,-100,100);
  drstcel->SetLineColor(kBlue);
  drshcel->SetLineColor(kRed);
  drstdel->SetLineColor(kBlue);
  drshdel->SetLineColor(kRed);
  drstceg->SetLineColor(kBlue);
  drshceg->SetLineColor(kRed);
  drstdeg->SetLineColor(kBlue);
  drshdeg->SetLineColor(kRed);

  spdiag->Project("drstcel","strho-mcshrho",stereo+ce+ldist);
  spdiag->Project("drstceg","strho-mcshrho",stereo+ce+gdist);
  spdiag->Project("drstdel","strho-mcshrho",stereo+de+ldist);
  spdiag->Project("drstdeg","strho-mcshrho",stereo+de+gdist);
  spdiag->Project("drshcel","shrho-mcshrho",stereo+ce+ldist);
  spdiag->Project("drshceg","shrho-mcshrho",stereo+ce+gdist);
  spdiag->Project("drshdel","shrho-mcshrho",stereo+de+ldist);
  spdiag->Project("drshdeg","shrho-mcshrho",stereo+de+gdist);

  TLegend* drleg = new TLegend(0.1,0.7,0.5,0.9);
  drleg->AddEntry(drstcel,"StereoHit Position","L");
  drleg->AddEntry(drshcel,"StrawHit Position","L");

  TCanvas* drcan = new TCanvas("drcan", "Delta rho can",900,900);
  drcan->Divide(2,2);
  drcan->cd(1);
  drstcel->Draw();
  drshcel->Draw("same");
  drleg->Draw();
  drcan->cd(2);
  drstdel->Draw();
  drshdel->Draw("same");
  drcan->cd(3);
  drstceg->Draw();
  drshceg->Draw("same");
  drcan->cd(4);
  drstdeg->Draw();
  drshdeg->Draw("same");

  TH1F* dpstcel = new TH1F("dpstcel","Reco - True azimuth, #Delta z< 40 mm, CE;#Delta #phi (mm)",100,-0.5,0.5);
  TH1F* dpstdel = new TH1F("dpstdel","Reco - True azimuth, #Delta z< 40 mm, #delta e^{-};#Delta #phi (mm)",100,-0.5,0.5);
  TH1F* dpshcel = new TH1F("dpshcel","Reco - True azimuth, #Delta z< 40 mm, CE;#Delta #phi (mm)",100,-0.5,0.5);
  TH1F* dpshdel = new TH1F("dpshdel","Reco - True azimuth, #Delta z< 40 mm, #delta e^{-};#Delta #phi (mm)",100,-0.5,0.5);
  TH1F* dpstceg = new TH1F("dpstceg","Reco - True azimuth, #Delta z> 40 mm, CE;#Delta #phi (mm)",100,-0.5,0.5);
  TH1F* dpstdeg = new TH1F("dpstdeg","Reco - True azimuth, #Delta z> 40 mm, #delta e^{-};#Delta #phi (mm)",100,-0.5,0.5);
  TH1F* dpshceg = new TH1F("dpshceg","Reco - True azimuth, #Delta z> 40 mm, CE;#Delta #phi (mm)",100,-0.5,0.5);
  TH1F* dpshdeg = new TH1F("dpshdeg","Reco - True azimuth, #Delta z> 40 mm, #delta e^{-};#Delta #phi (mm)",100,-0.5,0.5);
  dpstcel->SetLineColor(kBlue);
  dpshcel->SetLineColor(kRed);
  dpstdel->SetLineColor(kBlue);
  dpshdel->SetLineColor(kRed);
  dpstceg->SetLineColor(kBlue);
  dpshceg->SetLineColor(kRed);
  dpstdeg->SetLineColor(kBlue);
  dpshdeg->SetLineColor(kRed);

  spdiag->Project("dpstcel","stphi-mcshphi",stereo+ce+ldist);
  spdiag->Project("dpstceg","stphi-mcshphi",stereo+ce+gdist);
  spdiag->Project("dpstdel","stphi-mcshphi",stereo+de+ldist);
  spdiag->Project("dpstdeg","stphi-mcshphi",stereo+de+gdist);
  spdiag->Project("dpshcel","shphi-mcshphi",stereo+ce+ldist);
  spdiag->Project("dpshceg","shphi-mcshphi",stereo+ce+gdist);
  spdiag->Project("dpshdel","shphi-mcshphi",stereo+de+ldist);
  spdiag->Project("dpshdeg","shphi-mcshphi",stereo+de+gdist);

  TLegend* dpleg = new TLegend(0.1,0.7,0.5,0.9);
  dpleg->AddEntry(dpstcel,"StereoHit Position","L");
  dpleg->AddEntry(dpshcel,"StrawHit Position","L");

  TCanvas* dpcan = new TCanvas("dpcan", "Delta phi can",900,900);
  dpcan->Divide(2,2);
  dpcan->cd(1);
  dpstcel->Draw();
  dpshcel->Draw("same");
  dpleg->Draw();
  dpcan->cd(2);
  dpstdel->Draw();
  dpshdel->Draw("same");
  dpcan->cd(3);
  dpstceg->Draw();
  dpshceg->Draw("same");
  dpcan->cd(4);
  dpstdeg->Draw();
  dpshdeg->Draw("same");
}
