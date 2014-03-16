#include <algorithm>
#include "TTree.h"
#include "TCut.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TLine.h"
#include "TStyle.h"
#include <iostream>
using namespace std;

void StrawHitTest (TTree* hits, char* page="bcan" ) {

  TString spage(page);
  gStyle->SetOptStat(0);
  TCut conv("mcpdg==11&&mcgen==2");
  TCut oele("abs(mcpdg)==11&&mcgen!=2");
  TCut dio("mcpdg==11&&mcgen==6");
  TCut delta("mcpdg==11&&mcgen<0&&mcproc==17");
  TCut pconv("mcpdg==11&&mcgen<0&&mcproc==11");
  TCut compt("mcpdg==11&&mcgen<0&&mcproc==12");
  TCut proton("mcpdg==2212");
  TCut photon("mcpdg==22");
  TCut neutron("mcpdg==2112");
  TCut bkge("abs(mcpdg)==11&&mcgen<0");
  TCut bkgo("abs(mcpdg)!=11&&mcpdg!=2212");
  TCut xtalk("xtalk");
  TCut direct("!xtalk");
  TCut opart=!(conv||oele||proton);

  TCut goodevt("mcmom>100&&mctd<1&&mctd>0.577");
  TCut goodpeak("abs(tpeak-mct0-25)<30");
  if(spage =="scan"){
    TCanvas* scan = new TCanvas("scan","simulation",1200,800);
    TH1F* pdgid = new TH1F("pdgid","PDG Id",100,-200,2300);
    pdgid->SetStats(0);
    TH1F* gide = new TH1F("gide","Generator code",21,-1.5,19.5);
    TH1F* gidp = new TH1F("gidp","Generator code",21,-1.5,19.5);
    TH1F* pide = new TH1F("pide","Process code",101,-1.5,99.5);
    TH1F* pidp = new TH1F("pidp","Process code",101,-1.5,99.5);
    gide->SetLineColor(kRed);
    gidp->SetLineColor(kBlue);
    pide->SetLineColor(kRed);
    pidp->SetLineColor(kBlue);

    TH1F* nuconv = new TH1F("nuconv","N unique trkids in StrawHit",6,-0.5,5.5);
    TH1F* nudio = new TH1F("nudio","N unique trkids in StrawHit",6,-0.5,5.5);
    TH1F* nudel = new TH1F("nudel","N unique trkids in StrawHit",6,-0.5,5.5);
    TH1F* nup = new TH1F("nup","N unique trkids in StrawHit",6,-0.5,5.5);
    nuconv->SetLineColor(kRed);
    nudio->SetLineColor(kGreen);
    nudel->SetLineColor(kCyan);
    nup->SetLineColor(kBlue);
    TH2F* eve = new TH2F("eve","StrawHit Edep vs #Sigma MC Edep;MeV;MeV",100,0,0.3,100,0,0.3);
    eve->SetStats(0);

    hits->Project("pdgid","mcpdg");
    hits->Project("gide","mcgen","mcpdg==11");
    hits->Project("gidp","mcgen",proton);
    hits->Project("pide","mcproc",delta);
    hits->Project("pidp","mcproc",proton);
    hits->Project("nuconv","mcnunique",conv);
    hits->Project("nudio","mcnunique",dio);
    hits->Project("nudel","mcnunique",delta);
    hits->Project("nup","mcnunique",proton);

    TLegend* leg = new TLegend(0.5,0.5,0.8,0.8);
    leg->AddEntry(gide,"Electrons","l");
    leg->AddEntry(gidp,"Protons","l");

    TLegend* leg3 = new TLegend(0.5,0.5,0.8,0.8);    
    leg3->AddEntry(nuconv,"Conv. Electrons","l");
    leg3->AddEntry(nudio,"DIO Electrons","l");
    leg3->AddEntry(nudel,"Delta Electrons","l");
    leg3->AddEntry(nup,"Protons","l");

    scan->Clear();
    scan->Divide(2,2);
    scan->cd(1);
    //gPad->SetLogy();
    gidp->Draw();
    gide->Draw("same");
    scan->cd(2);
    //gPad->SetLogy();
    pide->Draw();
    pidp->Draw("same");
    leg->Draw();
    scan->cd(3);
    nup->Draw();
    nuconv->Draw("same");
    nudio->Draw("same");
    nudel->Draw("same");
    leg3->Draw();

    scan->cd(4);
    hits->Draw("edep:mcedep>>eve","","",10000);

  } else if(spage =="tcan"){
    THStack* tstack = new THStack("tc","Reco Hit Time by source;Hit Time (ns);Hits/event/ns");
    TH1F* ctime = new TH1F("ctime","Conversion Reco Hit Time",150,250,1750);
    TH1F* ptime = new TH1F("ptime","Proton Reco Hit Time",150,250,1750);
    TH1F* etime = new TH1F("etime","Electron Reco Hit Time",150,250,1750);
    TH1F* otime = new TH1F("otime","Other Reco Hit Time",150,250,1750);
    ctime->SetFillColor(kRed);
    ptime->SetFillColor(kBlack);
    etime->SetFillColor(kBlue);
    otime->SetFillColor(kGreen);

    hits->Project("otime","time",opart);
    otime->Scale(1e-4);
    tstack->Add(otime);
    hits->Project("ctime","time",conv);
    ctime->Scale(1e-4);
    tstack->Add(ctime);
    hits->Project("ptime","time",proton);
    ptime->Scale(1e-4);
    tstack->Add(ptime);
    hits->Project("etime","time",oele);
    etime->Scale(1e-4);
    tstack->Add(etime);
    
    cout << "other integral = " << otime->Integral() 
    << "CE inegral = " << ctime->Integral()
    << "P inegral = " << ptime->Integral()
    << "e inegral = " << etime->Integral() << endl;

    TLegend* tleg = new TLegend(.6,.7,.9,.9);
    char title[50];
    snprintf(title,50,"CE, #int=%4.0f",ctime->Integral()*10.0);
    tleg->AddEntry(ctime,title,"F");
    snprintf(title,50,"Proton, #int=%4.0f",ptime->Integral()*10.0);
    tleg->AddEntry(ptime,title,"F");
    snprintf(title,50,"Other e^{-}, #int=%4.0f",etime->Integral()*10.0);
    tleg->AddEntry(etime,title,"F");
    snprintf(title,50,"Other Particles, #int=%4.0f",otime->Integral()*10.0);
    tleg->AddEntry(otime,title,"F");

    TCanvas* tcan = new TCanvas("tcan","tcan",800,800);
    tstack->Draw();
    tleg->Draw();

  } else if(spage == "bcan"){

    TCanvas* bcan = new TCanvas("bcan","background",1200,800);

    TH1F* econv = new TH1F("econv","StrawHit EDep;MeV",200,-0.001,0.015);
    TH1F* edio = new TH1F("edio","StrawHit EDep;MeV",200,-0.001,0.015);
//    TH1F* eneut = new TH1F("eneut","StrawHit EDep;MeV",200,-0.001,0.015);
//    TH1F* ephot = new TH1F("ephot","StrawHit EDep;MeV",200,-0.001,0.015);
    TH1F* edelta = new TH1F("edelta","StrawHit EDep;MeV",200,-0.001,0.015);
    TH1F* ep = new TH1F("ep","StrawHit EDep;MeV",200,-0.001,0.015);
    TH1F* ex = new TH1F("ex","StrawHit EDep;MeV",200,-0.001,0.015);
    econv->SetLineColor(kRed);
    edio->SetLineColor(kGreen);
    edelta->SetLineColor(kCyan);
    ep->SetLineColor(kBlue);
    ex->SetLineColor(kMagenta);
//    eneut->SetLineColor(kYellow);
//    ephot->SetLineColor(kMagenta);
    econv->SetStats(0);
    edio->SetStats(0);
    edelta->SetStats(0);
    ep->SetStats(0);
    ex->SetStats(0);

    TH1F* rconv = new TH1F("rconv","StrawHit Radius;mm",100,360,750);
    TH1F* rdio = new TH1F("rdio","StrawHit Radius;mm",100,360,750);
//    TH1F* rneut = new TH1F("rneut","StrawHit Radius;mm",100,360,750);
//    TH1F* rphot = new TH1F("rphot","StrawHit Radius;mm",100,360,750);
    TH1F* rdelta = new TH1F("rdelta","StrawHit Radius;mm",100,360,750);
    TH1F* rp = new TH1F("rp","StrawHit Radius;mm",100,360,750);
    TH1F* rx = new TH1F("rx","StrawHit Radius;mm",100,360,750);
    rconv->SetLineColor(kRed);
    rp->SetMinimum(10);
    rdio->SetLineColor(kGreen);
    rdelta->SetLineColor(kCyan);
    rp->SetLineColor(kBlue);
    rx->SetLineColor(kMagenta);
//    rneut->SetLineColor(kYellow);
//    rphot->SetLineColor(kMagenta);
    /*   
	 TH1F* nconv = new TH1F("nconv","N, D<10cm",41,-1.5,39.5);
	 TH1F* ndio = new TH1F("ndio","N, D<10cm",41,-1.5,39.5);
	 TH1F* ndelta = new TH1F("ndelta","N, D<10cm",41,-1.5,39.5);
	 TH1F* np = new TH1F("np","N, D<10cm",41,-1.5,39.5);
	 nconv->SetLineColor(kRed);
	 ndio->SetLineColor(kGreen);
	 ndelta->SetLineColor(kCyan);
	 np->SetLineColor(kBlue);

	 TH1F* zconv = new TH1F("zconv","Z",100,-1500,1500);
	 TH1F* zdio = new TH1F("zdio","Z",100,-1500,1500);
	 TH1F* zdelta = new TH1F("zdelta","Z",100,-1500,1500);
	 TH1F* zp = new TH1F("zp","Z",100,-1500,1500);
	 zconv->SetLineColor(kRed);
	 zdio->SetLineColor(kGreen);
	 zdelta->SetLineColor(kCyan);
	 zp->SetLineColor(kBlue);
     */
    TCut ecut("edep<0.0045");
    TCut rmin("sqrt(shpos.x^2+shpos.y^2)>410");

    hits->Project("econv","edep",conv+direct+goodevt);
    hits->Project("edio","edep",dio+direct+goodevt);
//    hits->Project("eneut","edep",neutron+goodevt);
//    hits->Project("ephot","edep",photon+goodevt);
    hits->Project("edelta","edep",bkge+direct+goodevt);
    hits->Project("ep","edep",proton+direct+goodevt);
    hits->Project("ex","edep",xtalk+goodevt);

    hits->Project("rconv","sqrt(shpos.y^2+shpos.x^2)",conv+direct+goodevt);
    hits->Project("rdio","sqrt(shpos.y^2+shpos.x^2)",dio+direct+goodevt);
//    hits->Project("rneut","sqrt(shpos.y^2+shpos.x^2)",neutron+goodevt);
//    hits->Project("rphot","sqrt(shpos.y^2+shpos.x^2)",photon+goodevt);
    hits->Project("rdelta","sqrt(shpos.y^2+shpos.x^2)",bkge+direct+goodevt);
    hits->Project("rp","sqrt(shpos.y^2+shpos.x^2)",proton+direct+goodevt);
    hits->Project("rx","sqrt(shpos.y^2+shpos.x^2)",xtalk+goodevt);
    /*    
	  hits->Project("nconv","n200",conv+goodevt);
	  hits->Project("ndio","n200",bkge+goodevt);
	  hits->Project("ndelta","n200",bkgo+goodevt);
	  hits->Project("np","n200",proton+goodevt);

	  hits->Project("zconv","shpos.z",conv+goodevt);
	  hits->Project("zdio","shpos.z",bkge+goodevt);
	  hits->Project("zdelta","shpos.z",bkgo+goodevt);
	  hits->Project("zp","shpos.z",proton+goodevt);
     */

    bcan->Clear();
    bcan->Divide(1,2);
    bcan->cd(1);
    gPad->SetLogy();
    /*    edio->Draw();
	  econv->Draw("same");
	  ep->Draw("same");
	  edelta->Draw("same");
	  leg2->Draw();

	  bcan->cd(2);
     */
    edelta->GetXaxis()->SetRangeUser(-0.005,0.04);
    edelta->SetMinimum(econv->GetMaximum()/100.0);
    double maxv = edelta->GetMaximum();
    maxv = max(maxv, ep->GetMaximum());
    maxv = max(maxv, ex->GetMaximum());
    edelta->SetMaximum(1.1*maxv);
    edelta->Draw();
    ep->Draw("same");
    econv->Draw("same");
    edio->Draw("same");
    ex->Draw("same");
//    eneut->Draw("same");
//    ephot->Draw("same");

    TLine* ecut_t = new TLine(0.003,0.0,0.003,edelta->GetMaximum());
    ecut_t->SetLineColor(kBlack);
    ecut_t->SetLineStyle(2);
    ecut_t->SetLineWidth(2);
    TLine* ecut_l = new TLine(0.0055,0.0,0.0055,edelta->GetMaximum());
    ecut_l->SetLineColor(kBlack);
    ecut_l->SetLineStyle(3);
    ecut_l->SetLineWidth(2);
    ecut_t->Draw();
//    ecut_l->Draw();

    TLegend* leg2 = new TLegend(0.55,0.6,0.9,0.9);
    leg2->AddEntry(rconv,"Conv. Electrons","l");
    leg2->AddEntry(rdelta,"Bkg Electrons","l");
    leg2->AddEntry(rdio,"DIO Electrons","l");
//    leg2->AddEntry(rneut,"Neutrons","l");
//    leg2->AddEntry(rphot,"Photons","l");
    leg2->AddEntry(rp,"Protons","l");
    leg2->AddEntry(rx,"X-talk","l");
    leg2->Draw();
    TLegend* leg3 = new TLegend(0.55,0.55,0.9,0.6);
    leg3->AddEntry(ecut_t,"Tight cut","l");
//    leg3->AddEntry(ecut_l,"Loose cut","l");
    leg3->Draw();

    int istart = ex->FindFixBin(0.0);
    int istop = ex->FindFixBin(0.003);
    double xtint = ex->Integral(istart,istop);
    double ceint = econv->Integral(istart,istop);
    double pint = ep->Integral();
    double dint1 = edelta->Integral(istart,istop);
    double dint2 = edelta->Integral(istop,edelta->GetNbinsX());
    std::cout << "cross-talk integral = " << xtint
    << " conversion integral = " << ceint 
    << " proton integral = " << pint
    << " delta integral = " << dint1 << " " << dint2 << std::endl;

    bcan->cd(2);
    gPad->SetLogy();

    rdelta->SetMinimum(rconv->GetMaximum()/100.);
    rdelta->Draw();
    rconv->Draw("same");
    rp->Draw("same");
    rdio->Draw("same");
    rx->Draw("same");
//    rneut->Draw("same");
 //   rphot->Draw("same");

    TLine* rmin_t = new TLine(370,0.0,370,rdelta->GetMaximum());
    rmin_t->SetLineColor(kBlack);
    rmin_t->SetLineStyle(2);
    rmin_t->SetLineWidth(2);
    TLine* rmin_l = new TLine(390,0.0,390,rdelta->GetMaximum());
    rmin_l->SetLineColor(kBlack);
    rmin_l->SetLineStyle(3);
    rmin_l->SetLineWidth(2);

    TLine* rmax_t = new TLine(650,0.0,650,rdelta->GetMaximum());
    rmax_t->SetLineColor(kBlack);
    rmax_t->SetLineStyle(2);
    rmax_t->SetLineWidth(2);
    TLine* rmax_l = new TLine(650,0.0,650,rdelta->GetMaximum());
    rmax_l->SetLineColor(kBlack);
    rmax_l->SetLineStyle(3);
    rmax_l->SetLineWidth(2);
    rmin_t->Draw();
//    rmin_l->Draw();
    rmax_t->Draw();
//    rmax_l->Draw();

    /*    
	  bcan->cd(4);
	  gPad->SetLogy();

	  zdio->Draw();
	  zconv->Draw("same");
	  zp->Draw("same");
	  zdelta->Draw("same");
     */  
  } else if (spage == "ccan") {

    TCanvas* ccan = new TCanvas("ccan","cleaned hits",1200,800);
    TCut clean = goodpeak+TCut("tight>0&&delta==0");
        
	  TH1F* gid = new TH1F("gid","Generator code",21,-1.5,19.5);
	  TH1F* gidc = new TH1F("gidc","Generator code",21,-1.5,19.5);

	  TH1F* rres = new TH1F("rres","StrawHit Radius resolution;mm",100,-200,200);
   
    TH1F* pres = new TH1F("pres","StrawHit #phi resolution;mm",100,-0.5,0.5);
        

	  gid->SetLineColor(kBlue);
	  gidc->SetLineColor(kRed);
	  hits->Project("gid","mcgen",clean);
	  hits->Project("gidc","mcgen",conv);

	  hits->Project("rres","sqrt(shpos.y^2+shpos.x^2)-sqrt(mcshpos.y^2+mcshpos.x^2)");
    
    hits->Project("pres","atan2(shpos.y,shpos.x)-atan2(mcshpos.y,mcshpos.x)");

    /*
       TLegend* leg3 = new TLegend(0.4,0.75,0.7,0.9);
       leg3->AddEntry(gidc,"Conv. Electron hits","l");
       leg3->AddEntry(gid,"Selected hits","l");
       leg3->Draw();
     */
    gStyle->SetOptStat(0);
    ccan->Clear();
    ccan->Divide(1,2);
    /*    ccan->cd(2);
	  gidc->Draw();
	  gid->Draw("same");
	  ccan->cd(3);
	  rres->Fit("gaus");
     */
    ccan->cd(2);
    pres->Fit("gaus");

    TH1F* tconv = new TH1F("tconv","time WRT peak;nsec",101,-80,80);
    TH1F* tdio = new TH1F("tdio","time WRT peak;nsec",101,-80,80);
    TH1F* tdelta = new TH1F("tdelta","time WRT peak;nsec",101,-80,80);
    TH1F* tp = new TH1F("tp","time WRT peak;nsec",101,-80,80);
    tconv->SetLineColor(kRed);
    tdio->SetLineColor(kGreen);
    tdelta->SetLineColor(kCyan);
    tp->SetLineColor(kBlue);

    hits->Project("tconv","time-tpeak",conv+goodevt+clean);
    hits->Project("tdio","time-tpeak",bkge+goodevt+clean);
    hits->Project("tdelta","time-tpeak",bkgo+goodevt+clean);
    hits->Project("tp","time-tpeak",proton+goodevt+clean);

    TLegend* leg4 = new TLegend(0.65,0.6,0.9,0.95);
    leg4->AddEntry(tconv,"Conv. Electrons","l");
    leg4->AddEntry(tdio,"Bkg Electrons","l");
    leg4->AddEntry(tdelta,"Photons","l");
    leg4->AddEntry(tp,"Protons","l");


    TLine* tmin = new TLine(-45,0.0,-45,tconv->GetMaximum());
    tmin->SetLineColor(kBlack);
    tmin->SetLineStyle(2);
    tmin->SetLineWidth(2);
    TLine* tmax = new TLine(45,0.0,45,tconv->GetMaximum());
    tmax->SetLineColor(kBlack);
    tmax->SetLineStyle(2);
    tmax->SetLineWidth(2);

    double maxt = 1.1*std::max(tconv->GetMaximum(),tdio->GetMaximum());
    tconv->SetMaximum(maxt);

    ccan->cd(1);
    tconv->Draw();
    tdio->Draw("same");
    tdelta->Draw("same");
    tp->Draw("same");
    leg4->Draw();
    tmin->Draw();
    tmax->Draw();

    // compute efficiency and purity
    Float_t ncon = gidc->GetEntries();
    Float_t nsel = gid->GetEntries();
    Float_t nsdelta = gid->GetBinContent(1);
    Float_t nscon = gid->GetBinContent(4);
    Float_t nsdio = gid->GetBinContent(8);
    Float_t nsprot = gid->GetBinContent(14);
    cout << "selection efficiency = " << nscon/ncon << endl;
    cout << "selection purity = " << nscon/nsel << endl;
    cout << "bkg fraction deltas = " << nsdelta/nsel << " dio " << nsdio/nsel << " proton " << nsprot/nsel << endl;
  }
}
