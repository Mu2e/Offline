#include "TH1F.h"
#include "TF1.h"
#include "TTree.h"
#include "TLegend.h"
#include "TPad.h"
#include "TPaveText.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TLine.h"
#include "TArrow.h"
#include "TCut.h"
#include "TBox.h"
#include "TMath.h"
#include "TProfile.h"
#include "TDirectory.h"
#include "Math/Math.h"
#include "THStack.h"



void Cosmics(TTree* cr, const char* page="prod") {
  TString tpage(page);

  if(tpage=="prod"){
    TH2F* pxyep = new TH2F("pxyep","Electron production position;x(mm); y(mm)",100,-800,800,100,-800,800);
    TH2F* pxyem = new TH2F("pxyem","Electron production position;x(mm); y(mm)",100,-800,800,100,-800,800);
  
    pxyep->SetMarkerStyle(4);
    pxyep->SetMarkerColor(kBlue);
    pxyep->SetStats(0);
    pxyem->SetMarkerStyle(4);
    pxyem->SetMarkerColor(kRed);
    pxyem->SetStats(0);
    
    cr->Project("pxyem","ey:ex","mcpdgid==11");
    cr->Project("pxyep","ey:ex","mcpdgid==-11");

    TH2F* pxzep = new TH2F("pxzep","Electron production position;z(mm); x(mm)",100,1000,3600,100,-800,800);
    TH2F* pxzem = new TH2F("pxzem","Electron production position;z(mm); x(mm)",100,1000,3600,100,-800,800);
  
    pxzep->SetMarkerStyle(4);
    pxzep->SetMarkerColor(kBlue);
    pxzep->SetStats(0);
    pxzem->SetMarkerStyle(4);
    pxzem->SetMarkerColor(kRed);
    pxzem->SetStats(0);
  

    cr->Project("pxzem","ex:ez","mcpdgid==11");
    cr->Project("pxzep","ex:ez","mcpdgid==-11");

    TH2F* ppxzep = new TH2F("ppxzep","Projected parent production position y=0;z(mm); x(mm)",100,0,4500,100,-2000,2000);
    TH2F* ppxzem = new TH2F("ppxzem","Projected parent production position y=0;z(mm); x(mm)",100,0,4500,100,-2000,2000);
    TH2F* ppxzmp = new TH2F("ppxzmp","Projected parent production position y=0;z(mm); x(mm)",100,0,4500,100,-2000,2000);
    TH2F* ppxzmm = new TH2F("ppxzmm","Projected parent production position y=0;z(mm); x(mm)",100,0,4500,100,-2000,2000);
   
    ppxzep->SetMarkerStyle(4);
    ppxzep->SetMarkerColor(kBlue);
    ppxzep->SetLineColor(kBlue);
    ppxzep->SetStats(0);
    ppxzem->SetMarkerStyle(4);
    ppxzem->SetMarkerColor(kRed);
    ppxzem->SetLineColor(kRed);
    ppxzem->SetStats(0);
    ppxzmp->SetMarkerStyle(5);
    ppxzmp->SetMarkerColor(kCyan);
    ppxzmp->SetLineColor(kCyan);
    ppxzmp->SetStats(0);
    ppxzmm->SetMarkerStyle(5);
    ppxzmm->SetMarkerColor(kOrange);
    ppxzmm->SetLineColor(kOrange);
    ppxzmm->SetStats(0);
  
    cr->Project("ppxzem","ppx:ppz","mcpdgid==11");
    cr->Project("ppxzep","ppx:ppz","mcpdgid==-11");
    cr->Project("ppxzmm","ppx:ppz","mcpdgid==13");
    cr->Project("ppxzmp","ppx:ppz","mcpdgid==-13");

    TH1F* mcmomem = new TH1F("mcmomem","Tracker Momentum, reflected particle;MeV/c",100,25,250);
    TH1F* mcmomep = new TH1F("mcmomep","Tracker Momentum, reflected particle;MeV/c",100,25,250);
    TH1F* mcmommm = new TH1F("mcmommm","Tracker Momentum, reflected particle;MeV/c",100,25,250);
    TH1F* mcmommp = new TH1F("mcmommp","Tracker Momentum, reflected particle;MeV/c",100,25,250);
    mcmomem->SetStats(0);
    mcmomem->SetLineColor(kRed);
    mcmomep->SetStats(0);
    mcmomep->SetLineColor(kBlue);
    mcmommm->SetStats(0);
    mcmommm->SetLineColor(kOrange);
    mcmommp->SetStats(0);
    mcmommp->SetLineColor(kCyan);

    cr->Project("mcmomem","umcmom","mcpdgid==11");
    cr->Project("mcmomep","umcmom","mcpdgid==-11");
    cr->Project("mcmommm","umcmom","mcpdgid==13");
    cr->Project("mcmommp","umcmom","mcpdgid==-13");

    TH1F* mcd0em = new TH1F("mcd0em","Track D_{0}, reflected particle;mm",100,-500,500);
    TH1F* mcd0ep = new TH1F("mcd0ep","Track D_{0}, reflected particle;mm",100,-500,500);
    TH1F* mcd0mm = new TH1F("mcd0mm","Track D_{0}, reflected particle;mm",100,-500,500);
    TH1F* mcd0mp = new TH1F("mcd0mp","Track D_{0}, reflected particle;mm",100,-500,500);
    mcd0em->SetStats(0);
    mcd0em->SetLineColor(kRed);
    mcd0ep->SetStats(0);
    mcd0ep->SetLineColor(kBlue);
    mcd0mm->SetStats(0);
    mcd0mm->SetLineColor(kOrange);
    mcd0mp->SetStats(0);
    mcd0mp->SetLineColor(kCyan);

    cr->Project("mcd0em","umcd0","mcpdgid==11");
    cr->Project("mcd0ep","umcd0","mcpdgid==-11");
    cr->Project("mcd0mm","umcd0","mcpdgid==13");
    cr->Project("mcd0mp","umcd0","mcpdgid==-13");

    TH1F* pmomem = new TH1F("pmomem","Parent momentum;MeV/c",100,0,50000);
    TH1F* pmomep = new TH1F("pmomep","Parent momentum;MeV/c",100,0,50000);
    TH1F* pmommm = new TH1F("pmommm","Parent momentum;MeV/c",100,0,50000);
    TH1F* pmommp = new TH1F("pmommp","Parent momentum;MeV/c",100,0,50000);
    pmomem->SetStats(0);
    pmomem->SetLineColor(kRed);
    pmomep->SetStats(0);
    pmomep->SetLineColor(kBlue);
    pmommm->SetStats(0);
    pmommm->SetLineColor(kOrange);
    pmommp->SetStats(0);
    pmommp->SetLineColor(kCyan);

    cr->Project("pmomem","sqrt(pmx^2+pmy^2+pmz^2)","mcpdgid==11");
    cr->Project("pmomep","sqrt(pmx^2+pmy^2+pmz^2)","mcpdgid==-11");
    cr->Project("pmommm","sqrt(pmx^2+pmy^2+pmz^2)","mcpdgid==13");
    cr->Project("pmommp","sqrt(pmx^2+pmy^2+pmz^2)","mcpdgid==-13");

    TH1F* pmomlowem = new TH1F("pmomlowem","Parent momentum;MeV/c",100,0,10000);
    TH1F* pmomlowep = new TH1F("pmomlowep","Parent momentum;MeV/c",100,0,10000);
    TH1F* pmomlowmm = new TH1F("pmomlowmm","Parent momentum;MeV/c",100,0,10000);
    TH1F* pmomlowmp = new TH1F("pmomlowmp","Parent momentum;MeV/c",100,0,10000);
    pmomlowem->SetStats(0);
    pmomlowem->SetLineColor(kRed);
    pmomlowep->SetStats(0);
    pmomlowep->SetLineColor(kBlue);
    pmomlowmm->SetStats(0);
    pmomlowmm->SetLineColor(kOrange);
    pmomlowmp->SetStats(0);
    pmomlowmp->SetLineColor(kCyan);

    cr->Project("pmomlowem","sqrt(pmx^2+pmy^2+pmz^2)","mcpdgid==11");
    cr->Project("pmomlowep","sqrt(pmx^2+pmy^2+pmz^2)","mcpdgid==-11");
    cr->Project("pmomlowmm","sqrt(pmx^2+pmy^2+pmz^2)","mcpdgid==13");
    cr->Project("pmomlowmp","sqrt(pmx^2+pmy^2+pmz^2)","mcpdgid==-13");

    TH1F* pctem = new TH1F("pctem","Projected parent cos(#theta)",100,-1,1);
    TH1F* pctep = new TH1F("pctep","Projected parent cos(#theta)",100,-1,1);
    TH1F* pctmm = new TH1F("pctmm","Projected parent cos(#theta)",100,-1,1);
    TH1F* pctmp = new TH1F("pctmp","Projected parent cos(#theta)",100,-1,1);
    pctem->SetStats(0);
    pctem->SetLineColor(kRed);
    pctep->SetStats(0);
    pctep->SetLineColor(kBlue);
    pctmm->SetStats(0);
    pctmm->SetLineColor(kOrange);
    pctmp->SetStats(0);
    pctmp->SetLineColor(kCyan);
    cr->Project("pctem","pmz/sqrt(pmx^2+pmy^2+pmz^2)","mcpdgid==11");
    cr->Project("pctep","pmz/sqrt(pmx^2+pmy^2+pmz^2)","mcpdgid==-11");
    cr->Project("pctmm","pmz/sqrt(pmx^2+pmy^2+pmz^2)","mcpdgid==13");
    cr->Project("pctmp","pmz/sqrt(pmx^2+pmy^2+pmz^2)","mcpdgid==-13");
    
    TLegend* leg = new TLegend(0.7,0.5,0.9,0.9);
    leg->AddEntry(ppxzem,"e^{-}","LP");
    leg->AddEntry(ppxzep,"e^{+}","LP");
    leg->AddEntry(ppxzmm,"#mu^{-}","LP");
    leg->AddEntry(ppxzmp,"#mu^{+}","LP");

    TCanvas* pcan = new TCanvas("pcan","Production",1000,800);

    pcan->Divide(2,2);
    unsigned ican(1);
    pcan->cd(ican++);
    pmommp->Draw();
    pmommm->Draw("same");
    pmomem->Draw("same");
    pmomep->Draw("same");
    leg->Draw();

    pcan->cd(ican++);
    pmomlowmp->Draw();
    pmomlowmm->Draw("same");
    pmomlowem->Draw("same");
    pmomlowep->Draw("same");

    pcan->cd(ican++);
    pctem->Draw();
    pctep->Draw("same");
    pctmp->Draw("same");
    pctmm->Draw("same");

    pcan->cd(ican++);
    ppxzem->Draw();
    ppxzep->Draw("same");
    ppxzmm->Draw("same");
    ppxzmp->Draw("same");

    TCanvas* rcan = new TCanvas("rcan","Reflected",800,800);

    rcan->Divide(2,2);
    ican = 1;
    rcan->cd(ican++);
    mcmomem->Draw();
    mcmomep->Draw("same");
    mcmommp->Draw("same");
    mcmommm->Draw("same");
    leg->Draw();

    rcan->cd(ican++);
    mcd0em->Draw();
    mcd0ep->Draw("same");
    mcd0mp->Draw("same");
    mcd0mm->Draw("same");
 
    rcan->cd(ican++);
    pxyem->Draw();
    pxyep->Draw("same");

    rcan->cd(ican++);
    pxzem->Draw();
    pxzep->Draw("same");

    
  } else if(tpage =="select") {
    TH2F* cmom = new TH2F("cmom","Downstream vs Upstream momentum;P_{d} (MeV);P_{u} (MeV)",100,60,200,100,60,200);
    TH2F* ctand = new TH2F("ctand","Downstream vs Upstream tan#lambda;tan(#lambda)_{d};tan(#lambda)_{u}",100,-1.25,-0.25,100,0.25,1.25);
    TH2F* cp0 = new TH2F("cp0","Downstream vs Upstream #phi_{0};radians;radians",100,-3.15,3.15,100,-3.15,3.15);
    TH2F* cd0 = new TH2F("cd0","Downstream vs Upstream d_{0};mm;mm",100,-400,400,100,-400,400);
    TH1F* dt0e = new TH1F("dt0e","Downstream - Upstream t_{0};nsec",100,50,120);
    TH1F* dt0m = new TH1F("dt0m","Downstream - Upstream t_{0};nsec",100,50,120);
    cmom->SetStats(0);
    ctand->SetStats(0);
    cd0->SetStats(0);
    cp0->SetStats(0);
    dt0e->SetStats(0);
    dt0m->SetStats(0);
    dt0e->SetLineColor(kRed);
    dt0m->SetLineColor(kCyan);
    
    cr->Project("cmom","dmom:umom");
    cr->Project("ctand","dtd:utd");
    cr->Project("cp0","dp0:up0");
    cr->Project("cd0","dd0:ud0");
    cr->Project("dt0e","dt0-ut0","abs(mcpdgid)==11");
    cr->Project("dt0m","dt0-ut0","abs(mcpdgid)==13");

    TCanvas* scan = new TCanvas("scan","Selection",1200,800);
    scan->Divide(3,2);
    scan->cd(1);
    cmom->Draw();
    scan->cd(2);
    ctand->Draw();
    scan->cd(3);
    cp0->Draw();
    scan->cd(4);
    cd0->Draw();
    scan->cd(5);
    dt0e->Draw();
    dt0m->Draw("same");


  }



}
