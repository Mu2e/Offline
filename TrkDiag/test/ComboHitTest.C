#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCut.h"
#include "TLegend.h"
#include "TPaveStats.h"
#include <string>

void ComboHitTest(TTree* CHD, const char* page="count"){
  string spage(page);
  TCut Ce("mcgen==2&&mcproc==56");
  TCut Bkg("mcproc!=56&&mcpdg==11");
  TCut multi("nch>1");
  if(spage == "count" ) {
    TH1F* nshc = new TH1F("nshc","N Straw Hits",10,0.5,10.5);
    TH1F* nchc = new TH1F("nchc","N Combo Hits",10,0.5,10.5);
    TH1F* nshb = new TH1F("nshb","N Straw Hits",10,0.5,10.5);
    TH1F* nchb = new TH1F("nchb","N Combo Hits",10,0.5,10.5);
    nshc->SetLineColor(kRed);
    nshb->SetLineColor(kBlue);
    nchc->SetLineColor(kRed);
    nchb->SetLineColor(kBlue);
    CHD->Project("nshc","nsh",Ce);
    CHD->Project("nchc","nch",Ce);
    CHD->Project("nshb","nsh",Bkg);
    CHD->Project("nchb","nch",Bkg);
    nshc->Scale(nshb->GetEntries()/nshc->GetEntries());
    nchc->Scale(nchb->GetEntries()/nchc->GetEntries());
    TCanvas* ccan = new TCanvas("ccan","ccan",600,400);
    TLegend* cleg = new TLegend(0.5,0.3,0.8,0.5);
    cleg->AddEntry(nshb,"Background","L");
    cleg->AddEntry(nshc,"Ce (scaled)","L");
    ccan->Divide(2,1);
    ccan->cd(1);
    nshb->Draw();
    nshc->Draw("samesh");
    gPad->Update();
    TPaveStats *st = (TPaveStats*)nshb->FindObject("stats");
    Coord_t dx = st->GetX2NDC()-st->GetX1NDC();
    st->SetX2NDC(st->GetX1NDC());
    st->SetX1NDC(st->GetX1NDC()-dx);
    st->Draw();
    cleg->Draw();
    ccan->cd(2);
    nchb->Draw();
    nchc->Draw("samesh");
    gPad->Update();
    st = (TPaveStats*)nchb->FindObject("stats");
    dx = st->GetX2NDC()-st->GetX1NDC();
    st->SetX2NDC(st->GetX1NDC());
    st->SetX1NDC(st->GetX1NDC()-dx);
    st->Draw();
    cleg->Draw();
  } else if(spage =="wres"){
    TH1F* cwresc = new TH1F("cwresc","Multi-Hit Wire Distance Res;WD_{reco}-WD_{MC} (mm)",100,-300,300);
    TH1F* cwresb = new TH1F("cwresb","Multi-Hit Wire Distance Res;WD_{reco}-WD_{MC} (mm)",100,-300,300);
    TH1F* swresc = new TH1F("swresc","Single-Hit Wire Distance Res;WD_{reco}-WD_{MC} (mm)",100,-300,300);
    TH1F* swresb = new TH1F("swresb","Single-Hit Wire Distance Res;WD_{reco}-WD_{MC} (mm)",100,-300,300);
    TH1F* swress = new TH1F("swress","Single-Hit Wire Distance Res;WD_{reco}-WD_{MC} (mm)",100,-300,300);
    cwresc->SetLineColor(kRed);
    cwresb->SetLineColor(kBlue);
    swresc->SetLineColor(kRed);
    swresb->SetLineColor(kBlue);
    swress->SetLineColor(kGreen);
    CHD->Project("cwresc","wdist-mcdist",Ce+multi);
    CHD->Project("cwresb","wdist-mcdist",Bkg+multi);
    CHD->Project("swresc","wdist+_dwire-mcdist",Ce+multi);
    CHD->Project("swresb","wdist+_dwire-mcdist",Bkg+multi);
    CHD->Project("swress","wdist+_dwire-mcdist",!multi);
    cwresc->Scale(cwresb->GetEntries()/cwresc->GetEntries());
    swresc->Scale(swresb->GetEntries()/swresc->GetEntries());
    swress->Scale(swresb->GetEntries()/swress->GetEntries());
    TLegend* wleg = new TLegend(0.1,0.6,0.4,0.9);
    wleg->AddEntry(cwresb,"Background","L");
    wleg->AddEntry(cwresc,"Ce (scaled)","L");
    TCanvas* wrcan = new TCanvas("wrcan","wrcan",600,400);
    wrcan->Divide(2,1);
    wrcan->cd(1);
    cwresc->Fit("gaus","","",-60,60);
    cwresb->Draw("samehs");
    wleg->Draw();
    wrcan->cd(2);
    swresc->Fit("gaus","","",-80,80);
    swresb->Draw("samehs");
    swress->Draw("samehs");
  } else if(spage =="wpull"){
    TH1F* cwpullc = new TH1F("cwpullc","Multi-Hit Wire Distance Pull;(WD_{reco}-WD_{MC})/err",100,-10,10);
    TH1F* cwpullb = new TH1F("cwpullb","Multi-Hit Wire Distance Pull;(WD_{reco}-WD_{MC})/err",100,-10,10);
    TH1F* swpullc = new TH1F("swpullc","Single-Hit Wire Distance Pull;(WD_{reco}-WD_{MC})/err",100,-10,10);
    TH1F* swpullb = new TH1F("swpullb","Single-Hit Wire Distance Pull;(WD_{reco}-WD_{MC})/err",100,-10,10);
    TH1F* swpulls = new TH1F("swpulls","Single-Hit Wire Distance Pull;(WD_{reco}-WD_{MC})/err",100,-10,10);
    cwpullc->SetLineColor(kRed);
    cwpullb->SetLineColor(kBlue);
    swpullc->SetLineColor(kRed);
    swpullb->SetLineColor(kBlue);
    swpulls->SetLineColor(kGreen);
    CHD->Project("cwpullc","(wdist-mcdist)/wres",Ce+multi);
    CHD->Project("cwpullb","(wdist-mcdist)/wres",Bkg+multi);
    CHD->Project("swpullc","(wdist+_dwire-mcdist)/_dwerr",Ce+multi);
    CHD->Project("swpullb","(wdist+_dwire-mcdist)/_dwerr",Bkg+multi);
    CHD->Project("swpulls","(wdist+_dwire-mcdist)/dwerr",!multi);
    cwpullc->Scale(cwpullb->GetEntries()/cwpullc->GetEntries());
    swpullc->Scale(swpullb->GetEntries()/swpullc->GetEntries());
    swpulls->Scale(swpullb->GetEntries()/swpulls->GetEntries());
    TLegend* wleg = new TLegend(0.1,0.6,0.4,0.9);
    wleg->AddEntry(cwpullb,"Background","L");
    wleg->AddEntry(cwpullc,"Ce (scaled)","L");
    TCanvas* wpcan = new TCanvas("wpcan","wpcan",600,400);
    wpcan->Divide(2,1);
    wpcan->cd(1);
    cwpullc->Fit("gaus","","",-2,2);
    cwpullb->Draw("samehs");
    wleg->Draw();
    wpcan->cd(2);
    swpullc->Fit("gaus","","",-2,2);
    swpullb->Draw("samehs");
    swpulls->Draw("samehs");

  } else if(spage == "res") {
    TH1F* rresc = new TH1F("rresc","Radial Resolution;R_{reco}-R_{MC} (mm)",100,-100.0,100.0);
    TH1F* fresc = new TH1F("fresc","#phi Resolution;R_{MC}#times(#phi_{reco}-#phi_{MC}) (mm)",100,-300.0,300.0);
    TH1F* rresb = new TH1F("rresb","Radial Resolution;R_{reco}-R_{MC} (mm)",100,-100.0,100.0);
    TH1F* fresb = new TH1F("fresb","#phi Resolution;R_{MC}#times(#phi_{reco}-#phi_{MC}) (mm)",100,-300.0,300.0);
    rresc->SetLineColor(kRed);
    rresb->SetLineColor(kBlue);
    fresc->SetLineColor(kRed);
    fresb->SetLineColor(kBlue);
    CHD->Project("rresc","sqrt(pos.fCoordinates.fY^2+pos.fCoordinates.fX^2)-sqrt(mcpos.fCoordinates.fY^2+mcpos.fCoordinates.fX^2)",Ce);
    CHD->Project("fresc","sqrt(mcpos.fCoordinates.fY^2+mcpos.fCoordinates.fX^2)*(atan2(pos.fCoordinates.fY,pos.fCoordinates.fX)-atan2(mcpos.fCoordinates.fY,mcpos.fCoordinates.fX))",Ce);
    CHD->Project("rresb","sqrt(pos.fCoordinates.fY^2+pos.fCoordinates.fX^2)-sqrt(mcpos.fCoordinates.fY^2+mcpos.fCoordinates.fX^2)",Bkg);
    CHD->Project("fresb","sqrt(mcpos.fCoordinates.fY^2+mcpos.fCoordinates.fX^2)*(atan2(pos.fCoordinates.fY,pos.fCoordinates.fX)-atan2(mcpos.fCoordinates.fY,mcpos.fCoordinates.fX))",Bkg);
    rresc->Scale(rresb->GetEntries()/rresc->GetEntries());
    fresc->Scale(fresb->GetEntries()/fresc->GetEntries());
    double rmax = max(rresc->GetMaximum(),rresb->GetMaximum());
    double fmax = max(fresc->GetMaximum(),fresb->GetMaximum());
    rresc->SetMaximum(1.1*rmax);
    rresb->SetMaximum(1.1*rmax);
    fresc->SetMaximum(1.1*fmax);
    fresb->SetMaximum(1.1*fmax);
    TLegend* resleg = new TLegend(0.1,0.6,0.4,0.9);
    resleg->AddEntry(rresb,"Background","L");
    resleg->AddEntry(rresc,"Ce (scaled)","L");
    TCanvas* rescan = new TCanvas("rescan","rescan",800,600);
    rescan->Divide(2,1);
    rescan->cd(1);
    rresb->Draw("h");
    rresc->Draw("samesh");
    gPad->Update();
    TPaveStats *st = (TPaveStats*)rresc->FindObject("stats");
    Coord_t dx = st->GetX2NDC()-st->GetX1NDC();
    st->SetX2NDC(st->GetX1NDC());
    st->SetX1NDC(st->GetX1NDC()-dx);
    st->Draw();
    resleg->Draw();
    rescan->cd(2);
    fresb->Draw("h");
    fresc->Draw("samesh");
    gPad->Update();
    st = (TPaveStats*)fresc->FindObject("stats");
    dx = st->GetX2NDC()-st->GetX1NDC();
    st->SetX2NDC(st->GetX1NDC());
    st->SetX1NDC(st->GetX1NDC()-dx);
    st->Draw();
  } else if (spage == "edep") {
    TH1F* pedep = new TH1F("pedep","ComboHit Energy Deposit;EDep (KeV)",100,0,4.0); 
    TH1F* eedep = new TH1F("eedep","ComboHit Energy Deposit;EDep (KeV)",100,0,4.0); 
    pedep->SetLineColor(kBlue);
    eedep->SetLineColor(kRed);
    pedep->SetStats(0);
    eedep->SetStats(0);
    CHD->Project("pedep","1000.0*edep","mcpdg==2212");
    CHD->Project("eedep","1000.0*edep","mcpdg==11");
    double np =pedep->GetEntries();
    double ne =eedep->GetEntries();
    pedep->Scale(1.0/np);
    eedep->Scale(1.0/ne);
    TCanvas* cedep = new TCanvas("cedep","cedep",600,600);
    pedep->Draw();
    eedep->Draw("same");
    TLegend* edleg = new TLegend(0.2,0.7,0.4,0.9);
    edleg->AddEntry(pedep,"protons","L");
    edleg->AddEntry(eedep,"electrons","L");
    edleg->Draw();
  }
}
