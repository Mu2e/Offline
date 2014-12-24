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



void Cosmics(TTree* cr, const char* page="parent") {
  TString tpage(page);

  if(tpage=="electron"){
    TH2F* pxyep = new TH2F("pxyep","Electron Production Position;x(mm); y(mm)",100,-800,800,100,-800,800);
    TH2F* pxyem = new TH2F("pxyem","Electron Production Position;x(mm); y(mm)",100,-800,800,100,-800,800);
  
    pxyep->SetMarkerStyle(4);
    pxyep->SetMarkerColor(kBlue);
    pxyep->SetStats(0);
    pxyem->SetMarkerStyle(4);
    pxyem->SetMarkerColor(kRed);
    pxyem->SetStats(0);
    
    cr->Project("pxyem","opos.y:opos.x","mc.pdg==11");
    cr->Project("pxyep","opos.y:opos.x","mc.pdg==-11");

    TH2F* pxzep = new TH2F("pxzep","Electron Production Position;z(mm); x(mm)",100,1000,3600,100,-800,800);
    TH2F* pxzem = new TH2F("pxzem","Electron Production Position;z(mm); x(mm)",100,1000,3600,100,-800,800);
  
    pxzep->SetMarkerStyle(4);
    pxzep->SetMarkerColor(kBlue);
    pxzep->SetStats(0);
    pxzem->SetMarkerStyle(4);
    pxzem->SetMarkerColor(kRed);
    pxzem->SetStats(0);
  
    cr->Project("pxzem","opos.x:opos.z","mc.pdg==11");
    cr->Project("pxzep","opos.x:opos.z","mc.pdg==-11");

    TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
    leg->AddEntry(pxyem,"e^{-}","P");
    leg->AddEntry(pxyep,"e^{+}","P");

    TCanvas* ecan = new TCanvas("ecan","Electrons",800,500);
    ecan->Divide(2,1);

    unsigned ican(1);
    ecan->cd(ican++);
    pxyem->Draw();
    pxyep->Draw("same");

    ecan->cd(ican++);
    pxzem->Draw();
    pxzep->Draw("same");
    leg->Draw();

  } else if(tpage=="parent") {

    TH1F* ppxep = new TH1F("ppxep","Projected Parent x Position, y=0;x(mm)",100,-2500,2500);
    TH1F* ppxem = new TH1F("ppxem","Projected Parent x Position, y=0;x(mm)",100,-2500,2500);
    TH1F* ppxmp = new TH1F("ppxmp","Projected Parent x Position, y=0;x(mm)",100,-2500,2500);
    TH1F* ppxmm = new TH1F("ppxmm","Projected Parent x Position, y=0;x(mm)",100,-2500,2500);

    TH1F* ppzep = new TH1F("ppzep","Projected Parent z Position, y=0;z(mm)",100,-1000,5000);
    TH1F* ppzem = new TH1F("ppzem","Projected Parent z Position, y=0;z(mm)",100,-1000,5000);
    TH1F* ppzmp = new TH1F("ppzmp","Projected Parent z Position, y=0;z(mm)",100,-1000,5000);
    TH1F* ppzmm = new TH1F("ppzmm","Projected Parent z Position, y=0;z(mm)",100,-1000,5000);

//    TH2F* ppxzep = new TH2F("ppxzep","Projected parent production position y=0;z(mm); x(mm)",50,-500,5000,50,-2500,2500);
//    TH2F* ppxzem = new TH2F("ppxzem","Projected parent production position y=0;z(mm); x(mm)",50,-500,5000,50,-2500,2500);
//    TH2F* ppxzmp = new TH2F("ppxzmp","Projected parent production position y=0;z(mm); x(mm)",50,-500,5000,50,-2500,2500);
//    TH2F* ppxzmm = new TH2F("ppxzmm","Projected parent production position y=0;z(mm); x(mm)",50,-500,5000,50,-2500,2500);
 
    TH1F* pmomem = new TH1F("pmomem","Parent Momentum;MeV/c",200,-1000,10000);
    TH1F* pmomep = new TH1F("pmomep","Parent Momentum;MeV/c",200,-1000,10000);
    TH1F* pmommm = new TH1F("pmommm","Parent Momentum;MeV/c",200,-1000,10000);
    TH1F* pmommp = new TH1F("pmommp","Parent Momentum;MeV/c",200,-1000,10000);

//    TH1F* pmomlowem = new TH1F("pmomlowem","Parent momentum;MeV/c",100,0,10000);
//    TH1F* pmomlowep = new TH1F("pmomlowep","Parent momentum;MeV/c",100,0,10000);
//    TH1F* pmomlowmm = new TH1F("pmomlowmm","Parent momentum;MeV/c",100,0,10000);
//    TH1F* pmomlowmp = new TH1F("pmomlowmp","Parent momentum;MeV/c",100,0,10000);

    TH1F* pctem = new TH1F("pctem","Parent cos(#theta)",100,-1,1);
    TH1F* pctep = new TH1F("pctep","Parent cos(#theta)",100,-1,1);
    TH1F* pctmm = new TH1F("pctmm","Parent cos(#theta)",100,-1,1);
    TH1F* pctmp = new TH1F("pctmp","Parent cos(#theta)",100,-1,1);

    ppxem->SetStats(0);
    ppxem->SetLineColor(kRed);
    ppxep->SetStats(0);
    ppxep->SetLineColor(kBlue);
    ppxmm->SetStats(0);
    ppxmm->SetLineColor(kOrange);
    ppxmp->SetStats(0);
    ppxmp->SetLineColor(kCyan);

    ppzem->SetStats(0);
    ppzem->SetLineColor(kRed);
    ppzep->SetStats(0);
    ppzep->SetLineColor(kBlue);
    ppzmm->SetStats(0);
    ppzmm->SetLineColor(kOrange);
    ppzmp->SetStats(0);
    ppzmp->SetLineColor(kCyan);

//    ppxzep->SetMarkerStyle(4);
//    ppxzep->SetMarkerColor(kBlue);
//    ppxzep->SetLineColor(kBlue);
//    ppxzep->SetStats(0);
//    ppxzem->SetMarkerStyle(4);
//    ppxzem->SetMarkerColor(kRed);
//    ppxzem->SetLineColor(kRed);
//    ppxzem->SetStats(0);
//    ppxzmp->SetMarkerStyle(5);
//    ppxzmp->SetMarkerColor(kCyan);
//    ppxzmp->SetLineColor(kCyan);
//    ppxzmp->SetStats(0);
//    ppxzmm->SetMarkerStyle(5);
//    ppxzmm->SetMarkerColor(kOrange);
//    ppxzmm->SetLineColor(kOrange);
//    ppxzmm->SetStats(0);

    pmomem->SetStats(0);
    pmomem->SetLineColor(kRed);
    pmomep->SetStats(0);
    pmomep->SetLineColor(kBlue);
    pmommm->SetStats(0);
    pmommm->SetLineColor(kOrange);
    pmommp->SetStats(0);
    pmommp->SetLineColor(kCyan);

//    pmomlowem->SetStats(0);
//    pmomlowem->SetLineColor(kRed);
//    pmomlowep->SetStats(0);
//    pmomlowep->SetLineColor(kBlue);
//    pmomlowmm->SetStats(0);
//    pmomlowmm->SetLineColor(kOrange);
//    pmomlowmp->SetStats(0);
//    pmomlowmp->SetLineColor(kCyan);

    pctem->SetStats(0);
    pctem->SetLineColor(kRed);
    pctep->SetStats(0);
    pctep->SetLineColor(kBlue);
    pctmm->SetStats(0);
    pctmm->SetLineColor(kOrange);
    pctmp->SetStats(0);
    pctmp->SetLineColor(kCyan);

    cr->Project("ppxem","pppos.x","mc.pdg==11");
    cr->Project("ppxep","pppos.x","mc.pdg==-11");
    cr->Project("ppxmm","pppos.x","mc.pdg==13");
    cr->Project("ppxmp","pppos.x","mc.pdg==-13");
    cr->Project("ppzem","pppos.z","mc.pdg==11");
    cr->Project("ppzep","pppos.z","mc.pdg==-11");
    cr->Project("ppzmm","pppos.z","mc.pdg==13");
    cr->Project("ppzmp","pppos.z","mc.pdg==-13");

//    cr->Project("ppxzem","pppos.x:pppos.z","mc.pdg==11");
//    cr->Project("ppxzep","pppos.x:pppos.z","mc.pdg==-11");
//    cr->Project("ppxzmm","pppos.x:pppos.z","mc.pdg==13");
//    cr->Project("ppxzmp","pppos.x:pppos.z","mc.pdg==-13");

//    cr->Project("pmomlowem","sqrt(pmom.x^2+pmom.y^2+pmom.z^2)","mc.pdg==11");
//    cr->Project("pmomlowep","sqrt(pmom.x^2+pmom.y^2+pmom.z^2)","mc.pdg==-11");
//    cr->Project("pmomlowmm","sqrt(pmom.x^2+pmom.y^2+pmom.z^2)","mc.pdg==13");
//    cr->Project("pmomlowmp","sqrt(pmom.x^2+pmom.y^2+pmom.z^2)","mc.pdg==-13");

    cr->Project("pmomem","sqrt(pmom.x^2+pmom.y^2+pmom.z^2)","mc.pdg==11");
    cr->Project("pmomep","sqrt(pmom.x^2+pmom.y^2+pmom.z^2)","mc.pdg==-11");
    cr->Project("pmommm","sqrt(pmom.x^2+pmom.y^2+pmom.z^2)","mc.pdg==13");
    cr->Project("pmommp","sqrt(pmom.x^2+pmom.y^2+pmom.z^2)","mc.pdg==-13");

    cr->Project("pctem","pmom.z/sqrt(pmom.x^2+pmom.y^2+pmom.z^2)","mc.pdg==11");
    cr->Project("pctep","pmom.z/sqrt(pmom.x^2+pmom.y^2+pmom.z^2)","mc.pdg==-11");
    cr->Project("pctmm","pmom.z/sqrt(pmom.x^2+pmom.y^2+pmom.z^2)","mc.pdg==13");
    cr->Project("pctmp","pmom.z/sqrt(pmom.x^2+pmom.y^2+pmom.z^2)","mc.pdg==-13");
    
    TLine* momcut = new TLine(1500,0,1500,pmommm->GetMaximum());
    momcut->SetLineStyle(3);
    momcut->SetLineWidth(3);
    TArrow* momcutdir = new TArrow(500,pmommm->GetMaximum(), 1500, pmommm->GetMaximum(),0.02,"|>");
    momcutdir->SetLineWidth(3);
    momcutdir->SetFillColor(kWhite);

    TLine* costcut = new TLine(0,0,0,pctmm->GetMaximum());
    costcut->SetLineStyle(3);
    costcut->SetLineWidth(3);
    TArrow* costcutdir = new TArrow(0.2,pctmm->GetMaximum(), 0.0, pctmm->GetMaximum(),0.02,"|>");
    costcutdir->SetLineWidth(3);
    costcutdir->SetFillColor(kWhite);

    TLine* x1cut = new TLine(-1500,0,-1500,ppxmm->GetMaximum());
    x1cut->SetLineStyle(3);
    x1cut->SetLineWidth(3);
    TArrow* x1cutdir = new TArrow(-2000,ppxmm->GetMaximum(), -1500, ppxmm->GetMaximum(),0.02,"|>");
    x1cutdir->SetLineWidth(3);
    x1cutdir->SetFillColor(kWhite);

    TLine* x2cut = new TLine(1500,0,1500,ppxmm->GetMaximum());
    x2cut->SetLineStyle(3);
    x2cut->SetLineWidth(3);
    TArrow* x2cutdir = new TArrow(2000,ppxmm->GetMaximum(), 1500, ppxmm->GetMaximum(),0.02,"|>");
    x2cutdir->SetLineWidth(3);
    x2cutdir->SetFillColor(kWhite);

    TLine* z1cut = new TLine(800,0,800,ppzmm->GetMaximum());
    z1cut->SetLineStyle(3);
    z1cut->SetLineWidth(3);
    TArrow* z1cutdir = new TArrow(300,ppzmm->GetMaximum(), 800, ppzmm->GetMaximum(),0.02,"|>");
    z1cutdir->SetLineWidth(3);
    z1cutdir->SetFillColor(kWhite);

    TLine* z2cut = new TLine(3200,0,3200,ppzmm->GetMaximum());
    z2cut->SetLineStyle(3);
    z2cut->SetLineWidth(3);
    TArrow* z2cutdir = new TArrow(3700,ppzmm->GetMaximum(), 3200, ppzmm->GetMaximum(),0.02,"|>");
    z2cutdir->SetLineWidth(3);
    z2cutdir->SetFillColor(kWhite);

    TLegend* leg = new TLegend(0.6,0.6,0.9,0.9);
    leg->AddEntry(ppxem,"e^{-}","LP");
    leg->AddEntry(ppxep,"e^{+}","LP");
    leg->AddEntry(ppxmm,"#mu^{-}","LP");
    leg->AddEntry(ppxmp,"#mu^{+}","LP");
    leg->AddEntry(momcutdir,"Generator Cut","L");

    TCanvas* pcan = new TCanvas("pcan","Parent",1000,800);

    pcan->Divide(2,2);
    unsigned ican(1);
    pcan->cd(ican++);
    pmommp->Draw();
    pmommm->Draw("same");
    pmomem->Draw("same");
    pmomep->Draw("same");
    leg->Draw();
    momcut->Draw();
    momcutdir->Draw();

//    pcan->cd(ican++);
//    pmomlowmp->Draw();
//    pmomlowmm->Draw("same");
//    pmomlowem->Draw("same");
//    pmomlowep->Draw("same");

    pcan->cd(ican++);
    pctem->Draw();
    pctep->Draw("same");
    pctmp->Draw("same");
    pctmm->Draw("same");
    costcut->Draw();
    costcutdir->Draw();


    pcan->cd(ican++);
    ppxem->Draw();
    ppxep->Draw("same");
    ppxmm->Draw("same");
    ppxmp->Draw("same");
    x1cut->Draw();
    x1cutdir->Draw();
    x2cut->Draw();
    x2cutdir->Draw();


    pcan->cd(ican++);
    ppzem->Draw();
    ppzep->Draw("same");
    ppzmm->Draw("same");
    ppzmp->Draw("same");
    z1cut->Draw();
    z1cutdir->Draw();
    z2cut->Draw();
    z2cutdir->Draw();

  } else if (tpage=="reflected") {
    TH1F* mcmomem = new TH1F("mcmomem","Reflected Particle Momentum;MeV/c",100,25,250);
    TH1F* mcmomep = new TH1F("mcmomep","Reflected Particle Momentum;MeV/c",100,25,250);
    TH1F* mcmommm = new TH1F("mcmommm","Reflected Particle Momentum;MeV/c",100,25,250);
    TH1F* mcmommp = new TH1F("mcmommp","Reflected Particle Momentum;MeV/c",100,25,250);

    TH1F* mcd0em = new TH1F("mcd0em","Reflected Particle d_{0};mm",100,-500,500);
    TH1F* mcd0ep = new TH1F("mcd0ep","Reflected Particle d_{0};mm",100,-500,500);
    TH1F* mcd0mm = new TH1F("mcd0mm","Reflected Particle d_{0};mm",100,-500,500);
    TH1F* mcd0mp = new TH1F("mcd0mp","Reflected Particle d_{0};mm",100,-500,500);
    mcmomem->SetStats(0);
    mcmomem->SetLineColor(kRed);
    mcmomep->SetStats(0);
    mcmomep->SetLineColor(kBlue);
    mcmommm->SetStats(0);
    mcmommm->SetLineColor(kOrange);
    mcmommp->SetStats(0);
    mcmommp->SetLineColor(kCyan);

    mcd0em->SetStats(0);
    mcd0em->SetLineColor(kRed);
    mcd0ep->SetStats(0);
    mcd0ep->SetLineColor(kBlue);
    mcd0mm->SetStats(0);
    mcd0mm->SetLineColor(kOrange);
    mcd0mp->SetStats(0);
    mcd0mp->SetLineColor(kCyan);

    cr->Project("mcmomem","umcent.mom","mc.pdg==11");
    cr->Project("mcmomep","umcent.mom","mc.pdg==-11");
    cr->Project("mcmommm","umcent.mom","mc.pdg==13");
    cr->Project("mcmommp","umcent.mom","mc.pdg==-13");

    cr->Project("mcd0em","umcent.d0","mc.pdg==11");
    cr->Project("mcd0ep","umcent.d0","mc.pdg==-11");
    cr->Project("mcd0mm","umcent.d0","mc.pdg==13");
    cr->Project("mcd0mp","umcent.d0","mc.pdg==-13");

    TLegend* leg = new TLegend(0.7,0.5,0.9,0.9);
    leg->AddEntry(mcmomem,"e^{-}","LP");
    leg->AddEntry(mcmomep,"e^{+}","LP");
    leg->AddEntry(mcmommm,"#mu^{-}","LP");
    leg->AddEntry(mcmommp,"#mu^{+}","LP");

    TCanvas* rcan = new TCanvas("rcan","Reflected",800,800);

    rcan->Divide(2,2);
    unsigned ican = 1;
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

   
  } else if(tpage =="pselect") {
    TH2F* cmom = new TH2F("cmom","Downstream vs Upstream momentum;P_{d} (MeV);P_{u} (MeV)",100,30,250,100,30,250);
    TH2F* ctand = new TH2F("ctand","Downstream vs Upstream tan#lambda;tan(#lambda)_{d};tan(#lambda)_{u}",100,-1.25,-0.25,100,0.25,1.25);
    TH2F* cp0 = new TH2F("cp0","Downstream vs Upstream #phi_{0};radians;radians",100,-3.15,3.15,100,-3.15,3.15);
    TH2F* cd0 = new TH2F("cd0","Downstream vs Upstream d_{0};mm;mm",100,-600,600,100,-600,600);
    cmom->SetStats(0);
    ctand->SetStats(0);
    cd0->SetStats(0);
    cp0->SetStats(0);
   
    cr->Project("cmom","dtrk.mom:utrk.mom");
    cr->Project("ctand","dtrk.td:utrk.td");
    cr->Project("cp0","dtrk.p0:utrk.p0");
    cr->Project("cd0","dtrk.d0:utrk.d0");
    TCanvas* scan = new TCanvas("scan","Selection",1200,800);
    scan->Divide(2,2);
    scan->cd(1);
    cmom->Draw("colorz");
    scan->cd(2);
    ctand->Draw("colorz");
    scan->cd(3);
    cp0->Draw("colorz");
    scan->cd(4);
    cd0->Draw("colorz");

  } else if(tpage == "dt") {
    TH1F* dtrkt0e = new TH1F("dtrkt0e","Downstream - Upstream t_{0};nsec",100,50,200);
    TH1F* dtrkt0m = new TH1F("dtrkt0m","Downstream - Upstream t_{0};nsec",100,50,200);
    dtrkt0e->SetLineColor(kRed);
    dtrkt0m->SetLineColor(kCyan);
    dtrkt0e->SetStats(0);
    dtrkt0m->SetStats(0);
    cr->Project("dtrkt0e","dtrk.t0-utrk.t0","abs(mc.pdg)==11");
    cr->Project("dtrkt0m","dtrk.t0-utrk.t0","abs(mc.pdg)==13");
    TLegend* leg = new TLegend(0.7,0.5,0.9,0.9);
    leg->AddEntry(dtrkt0e,"e^{#pm}","L");
    leg->AddEntry(dtrkt0m,"#mu^{#pm}","L");
    TCanvas* dt = new TCanvas("dt","Time difference",600,600);
    dtrkt0e->Draw();
    dtrkt0m->Draw("same");
    leg->Draw();

    TLine* ecutl = new TLine(60,0,60,dtrkt0e->GetMaximum());
    ecutl->SetLineStyle(3);
    ecutl->SetLineWidth(3);
    TArrow* ecutldir = new TArrow(50,dtrkt0e->GetMaximum(), 60, dtrkt0e->GetMaximum(),0.02,"|>");
    ecutldir->SetLineWidth(3);
    ecutldir->SetFillColor(kRed);
    ecutl->Draw();
    ecutldir->Draw();

    TLine* ecuth = new TLine(90,0,90,dtrkt0e->GetMaximum());
    ecuth->SetLineStyle(3);
    ecuth->SetLineWidth(3);
    TArrow* ecuthdir = new TArrow(100,dtrkt0e->GetMaximum(), 90, dtrkt0e->GetMaximum(),0.02,"|>");
    ecuthdir->SetLineWidth(3);
    ecuthdir->SetFillColor(kRed);
    ecuth->Draw();
    ecuthdir->Draw();

    TLine* mucutl = new TLine(75,0,75,dtrkt0m->GetMaximum());
    mucutl->SetLineStyle(3);
    mucutl->SetLineWidth(3);
    TArrow* mucutldir = new TArrow(65,dtrkt0m->GetMaximum(), 75, dtrkt0m->GetMaximum(),0.02,"|>");
    mucutldir->SetLineWidth(3);
    mucutldir->SetFillColor(kCyan);
    mucutl->Draw();
    mucutldir->Draw();

    TLine* mucuth = new TLine(150,0,150,dtrkt0m->GetMaximum());
    mucuth->SetLineStyle(3);
    mucuth->SetLineWidth(3);
    TArrow* mucuthdir = new TArrow(160,dtrkt0m->GetMaximum(), 150, dtrkt0m->GetMaximum(),0.02,"|>");
    mucuthdir->SetLineWidth(3);
    mucuthdir->SetFillColor(kCyan);
    mucuth->Draw();
    mucuthdir->Draw();


    
  } else if(tpage == "momdiff" ){

    TH2F* umomdiff = new TH2F("umomdiff","Upstream particle #Delta P, tracker exit - entrance;True #Delta P (MeV/c);Reco #Delta P (MeV/c)",100,-1,5,100,-1,5);
    TH2F* dmomdiff = new TH2F("dmomdiff","Downstream particle #Delta P, tracker entrance - exit;True #Delta P (MeV/c);Reco #Delta P (MeV/c)",100,-1,5,100,-1,5);
    umomdiff->SetStats(0);
    dmomdiff->SetStats(0);
    cr->Project("umomdiff","utrkxit.mom-utrk.mom:umcxit.mom-umcent.mom","umcent.mom>100&&utrk.pdg==mc.pdg");
    cr->Project("dmomdiff","dtrk.mom-dtrkxit.mom:dmcent.mom-dmcxit.mom","dmcent.mom>100&&dtrk.pdg==mc.pdg");
    TCanvas* mdcan = new TCanvas("mdcan","Momentum Difference",1000,800);
    TLine* diag = new TLine(0.0,0.0,5.0,5.0);
    diag->SetLineColor(kBlack);
    diag->SetLineStyle(2);
    diag->SetLineWidth(2);
    mdcan->Divide(2,1);
    mdcan->cd(1);
    umomdiff->Draw("colorz");
    diag->Draw();
    mdcan->cd(2);
    dmomdiff->Draw("colorz");
    diag->Draw();

  }
}
