#include <algorithm>
#include "TTree.h"
#include "TCut.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TLine.h"
#include "TBox.h"
#include "TColor.h"
#include "TStyle.h"
#include <iostream>
using namespace std;

// compactify codes
Int_t myorigin(Int_t part, Int_t gen) {
  enum mp{ beamproton=0, proton=1, photon=2, neutron=3, electron=4, ootmuon=5, other=6};
  if(part==2212&&gen==16)return beamproton;
  if(part==2212&&gen==28)return proton;
  if(part==22) return photon;
  if(part==11)return electron;
  if(part==2112)return neutron;
  if(part==13)return ootmuon;
  return other;
}


Int_t mypart(Int_t part) {
  enum mp{ proton=0, electron=1, photon=2, neutron=3, positron=4, muon=5, other=6};
  if(part==2212)return proton;
  if(part==11)return electron;
  if(part==22) return photon;
  if(part==-11)return positron;
  if(part==2112)return neutron;
  if(part==13 || part==-13)return muon;
  return other;
}

Int_t myppart(Int_t ppart, Int_t part){
  if(ppart==0)return mypart(part);
  return mypart(ppart);
}

Int_t myproc(Int_t proc,Bool_t xtalk) {
  enum mp{ primary=0, compt=1, delta=2, gamconv=3, photo=4, had=5, xt=6, other=7};
  if(xtalk)return xt;
  if(proc==56)return primary;
  if(proc==12)return compt;
  if(proc==17 || proc==21) return delta;
  if(proc==13)return gamconv;
  if(proc==40)return photo;
  if(proc==58)return had;
  return other;
}

Int_t myhpart(Int_t mcpdg,Int_t mcgen, Int_t mcproc){
  enum hpart{Proton=0,LowEe,DIO,Other,CE};
  if(mcpdg==11&&mcgen==2)return CE;
  if(mcpdg==11&&mcproc==56)return DIO;
  if(mcpdg==2212)return Proton;
  if(mcpdg==11)return LowEe;
  return Other;
}

void StrawHitTest (TTree* hits, const char* page="bcan",unsigned nevents=1000 ) {

  TString spage(page);
  gStyle->SetOptStat(0);
//  TCut conv("mcpdg==11&&mcgen==2&&mcmom>100.0");
  TCut conv("mcpdg==11&&mcgen==2");
  TCut oele("abs(mcpdg)==11&&mcgen!=2");
  TCut dio("mcpdg==11&&mcgen==6");
  TCut delta("mcpdg==11&&mcgen<0&&mcproc==17");
  TCut pconv("mcpdg==11&&mcgen<0&&mcproc==11");
  TCut compt("mcpdg==11&&mcgen<0&&mcproc==12");
  TCut proton("mcpdg==2212");
  TCut photon("mcpdg==22");
  TCut neutron("mcpdg==2112");
  TCut bkge("abs(mcpdg)==11&&mcgen!=2");
  TCut bkgo("abs(mcpdg)!=11&&mcpdg!=2212");
  TCut xtalk("mcxtalk");
  TCut direct("!mcxtalk");
  TCut opart=!(conv||oele||proton);

  TCut pproton("mcppdg==2212");
  TCut pneuton("mcppdg==2112");
  TCut pele("mcppdg==11");
  TCut ppos("mcppdg===11");
  TCut pgam("mcppdg==22");
  TCut nopar("mcppdg==0");

  TCut dioorigin("mcgpdg==11&&mcgid==28");
  TCut ootmuonorigin("mcgpdg==13");
  TCut norigin("mcgpdg==2112");
  TCut porigin("mcgpdg==22");
  TCut stpprotonorigin("mcgpdg==2212&&mcgid==28");
  TCut pprotonorigin("mcgpdg==2212&&mcgid==16");

  TCut hitsel("esel&&rsel&&tsel&&(!delta)&&(!isolated)");

  TCut goodevt("mcom>100");
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
    TH2F* eve = new TH2F("eve","StrawHit Edep vs #Sigma MC Edep;MeV;MeV",100,0,0.05,100,0,0.05);
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

  } else if(spage=="particle"){

    double pscale(1.0/nevents);
    TH1F* part = new TH1F("part","Particle Producing Hits;Particle;Hits/event",7,-0.5,6.5);
    hits->Project("part","mypart(mcpdg)");
    part->Scale(pscale);

    TH1F* ppart = new TH1F("ppart","Immediate Parent Particle Producing Hits;Parent Particle;Hits/event",7,-0.5,6.5);
    hits->Project("ppart","myppart(mcppdg,mcpdg)");
    ppart->Scale(pscale);

    TH1F* upart = new TH1F("upart","Ultimate Parent Particle Producing Hits;Ultimate Parent;Hits/event",7,-0.5,6.5);
    hits->Project("upart","myorigin(mcgpdg,mcgid)");
    upart->Scale(pscale);
 
    TH1F* proc = new TH1F("proc","Production Process of Particle Producing Hits;Production Process;Hits/event",8,-0.5,7.5);
    hits->Project("proc","myproc(mcproc,xtalk)");
    proc->Scale(pscale);
    
    TAxis* xax(0);
    int ibin(1);

    xax = part->GetXaxis();
    xax->SetBinLabel(ibin++,"Proton");
    xax->SetBinLabel(ibin++,"Electron");
    xax->SetBinLabel(ibin++,"Photon");
    xax->SetBinLabel(ibin++,"Neutron");
    xax->SetBinLabel(ibin++,"Positron");
    xax->SetBinLabel(ibin++,"Muon");
    xax->SetBinLabel(ibin++,"Other");

    ibin=1;
    xax = ppart->GetXaxis();
    xax->SetBinLabel(ibin++,"Proton");
    xax->SetBinLabel(ibin++,"Electron");
    xax->SetBinLabel(ibin++,"Photon");
    xax->SetBinLabel(ibin++,"Neutron");
    xax->SetBinLabel(ibin++,"Positron");
    xax->SetBinLabel(ibin++,"Muon");
    xax->SetBinLabel(ibin++,"Other");

    ibin=1;
    xax = upart->GetXaxis();
    xax->SetBinLabel(ibin++,"Beam Proton");
    xax->SetBinLabel(ibin++,"Target Proton");
    xax->SetBinLabel(ibin++,"Target Photon");
    xax->SetBinLabel(ibin++,"Target Neutron");
    xax->SetBinLabel(ibin++,"Target Electron");
    xax->SetBinLabel(ibin++,"Out-Of-Target");
    xax->SetBinLabel(ibin++,"Other");

    ibin=1;
    xax = proc->GetXaxis();
    xax->SetBinLabel(ibin++,"Primary");
    xax->SetBinLabel(ibin++,"Compton");
    xax->SetBinLabel(ibin++,"Delta-Ray");
    xax->SetBinLabel(ibin++,"#gamma Conversion");
    xax->SetBinLabel(ibin++,"Photoelectric");
    xax->SetBinLabel(ibin++,"Hadronic");
    xax->SetBinLabel(ibin++,"Cross-Talk");
    xax->SetBinLabel(ibin++,"Other");

    TCanvas* pcan = new TCanvas("pcan","pcan",1200,800);
    pcan->Divide(2,2);
    pcan->cd(1);
    part->Draw();
    pcan->cd(2);
    ppart->Draw();
    pcan->cd(3);
    upart->Draw();
    pcan->cd(4);
    proc->Draw();

  } else if(spage =="tcan"){
    THStack* tstack = new THStack("tc","Reco Hit Time by Particle;Hit Time (ns);Hits/event/ns");
    TH1F* ctime = new TH1F("ctime","Conversion Reco Hit Time",150,250,1750);
    TH1F* ptime = new TH1F("ptime","Proton Reco Hit Time",150,250,1750);
    TH1F* etime = new TH1F("etime","Electron Reco Hit Time",150,250,1750);
    TH1F* otime = new TH1F("otime","Other Reco Hit Time",150,250,1750);
    ctime->SetFillColor(kRed);
    ptime->SetFillColor(kBlack);
    etime->SetFillColor(kBlue);
    otime->SetFillColor(kGreen);

    double scale = 0.1/nevents;
    hits->Project("otime","time",opart);
    otime->Scale(scale);
    tstack->Add(otime);
    hits->Project("ctime","time",conv);
    ctime->Scale(scale);
    tstack->Add(ctime);
    hits->Project("ptime","time",proton);
    ptime->Scale(scale);
    tstack->Add(ptime);
    hits->Project("etime","time",oele);
    etime->Scale(scale);
    tstack->Add(etime);

    double total = otime->Integral() + ctime->Integral() + ptime->Integral() + etime->Integral();
 
    cout << "All other integral = " << otime->Integral() 
    << "CE inegral = " << ctime->Integral()
    << "P inegral = " << ptime->Integral()
    << "e inegral = " << etime->Integral() << " total = " << total << endl;

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
    snprintf(title,50,"Total #int=%4.0f",total*10.0);
    tleg->AddEntry((TObject*)0,title,"");


    TCanvas* tcan = new TCanvas("tcan","tcan",800,600);
    tstack->Draw("h");
    tleg->Draw();

  } else if(spage=="selcan"){

    THStack* tstacks = new THStack("tcs","Selected Reco Hit Time by Particle;Hit Time (ns);Hits/event/ns");
    TH1F* ctimes = new TH1F("ctimes","Conversion Reco Hit Time",150,250,1750);
    TH1F* ptimes = new TH1F("ptimes","Proton Reco Hit Time",150,250,1750);
    TH1F* etimes = new TH1F("etimes","Electron Reco Hit Time",150,250,1750);
    TH1F* otimes = new TH1F("otimes","Other Reco Hit Time",150,250,1750);
    ctimes->SetFillColor(kRed);
    ptimes->SetFillColor(kBlack);
    etimes->SetFillColor(kBlue);
    otimes->SetFillColor(kGreen);

    double scale = 0.1/nevents;
    hits->Project("otimes","time",hitsel+opart);
    otimes->Scale(scale);
    tstacks->Add(otimes);
    hits->Project("ctimes","time",hitsel+conv);
    ctimes->Scale(scale);
    tstacks->Add(ctimes);
    hits->Project("ptimes","time",hitsel+proton);
    ptimes->Scale(scale);
    tstacks->Add(ptimes);
    hits->Project("etimes","time",hitsel+oele);
    etimes->Scale(scale);
    tstacks->Add(etimes);
    
    cout << "Selected other integral = " << otimes->Integral() 
    << "CE inegral = " << ctimes->Integral()
    << "P inegral = " << ptimes->Integral()
    << "e inegral = " << etimes->Integral() << endl;

    TLegend* tlegs = new TLegend(.6,.7,.9,.9);
    char title[50];
    snprintf(title,50,"CE, #int=%4.0f",ctimes->Integral()*10.0);
    tlegs->AddEntry(ctimes,title,"F");
    snprintf(title,50,"Proton, #int=%4.0f",ptimes->Integral()*10.0);
    tlegs->AddEntry(ptimes,title,"F");
    snprintf(title,50,"Other e^{-}, #int=%4.0f",etimes->Integral()*10.0);
    tlegs->AddEntry(etimes,title,"F");
    snprintf(title,50,"Other Particles, #int=%4.0f",otimes->Integral()*10.0);
    tlegs->AddEntry(otimes,title,"F");

    TCanvas* tscan = new TCanvas("tscan","tscan",800,600);
    tscan->Divide(2,1);
    tscan->cd(1);
    tstacks->Draw("h");
    tlegs->Draw();

  } else if(spage == "bcan"){

    TH1F* econv = new TH1F("econv","Straw Hit Energy;Deposited Energy (KeV)",200,-1.0,16.0);
//    TH1F* edio = new TH1F("edio","Straw Hit Energy;Deposited Energy (KeV)",200,-1.0,16.0);
//    TH1F* eneut = new TH1F("eneut","Straw Hit Energy;Deposited Energy (KeV)",200,-1.0,16.0);
//    TH1F* ephot = new TH1F("ephot","Straw Hit Energy;Deposited Energy (KeV)",200,-1.0,16.0);
    TH1F* edelta = new TH1F("edelta","Straw Hit Energy;Deposited Energy (KeV)",200,-1.0,16.0);
    TH1F* ep = new TH1F("ep","Straw Hit Energy;Deposited Energy (KeV)",200,-1.0,16.0);
    TH1F* ex = new TH1F("ex","Straw Hit Energy;Deposited Energy (KeV)",200,-1.0,16.0);
    econv->SetLineColor(kRed);
//    edio->SetLineColor(kGreen);
    edelta->SetLineColor(kCyan);
    ep->SetLineColor(kBlue);
    ex->SetLineColor(kMagenta);
//    eneut->SetLineColor(kYellow);
//    ephot->SetLineColor(kMagenta);
    econv->SetStats(0);
//    edio->SetStats(0);
    edelta->SetStats(0);
    ep->SetStats(0);
    ex->SetStats(0);

    TH1F* rconv = new TH1F("rconv","Straw Hit Radius;Transverse Radius (mm)",100,360,720);
//    TH1F* rdio = new TH1F("rdio","Straw Hit Radius;Transverse Radius (mm)",100,360,720);
//    TH1F* rneut = new TH1F("rneut","Straw Hit Radius;Transverse Radius (mm)",100,360,720);
//    TH1F* rphot = new TH1F("rphot","Straw Hit Radius;Transverse Radius (mm)",100,360,720);
    TH1F* rdelta = new TH1F("rdelta","Straw Hit Radius;Transverse Radius (mm)",100,360,720);
    TH1F* rp = new TH1F("rp","Straw Hit Radius;Transverse Radius (mm)",100,360,720);
//    TH1F* rx = new TH1F("rx","Straw Hit Radius;Transverse Radius (mm)",100,360,720);
    rconv->SetLineColor(kRed);
    rp->SetMinimum(1);
//    rdio->SetLineColor(kGreen);
    rdelta->SetLineColor(kCyan);
    rp->SetLineColor(kBlue);
//    rx->SetLineColor(kMagenta);
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
    hits->Project("econv","edep*1000.0",conv+direct);
//    hits->Project("edio","edep*1000.0",dio+direct);
//    hits->Project("eneut","edep*1000.0",neutron);
//    hits->Project("ephot","edep*1000.0",photon);
    hits->Project("edelta","edep*1000.0",bkge+direct);
    hits->Project("ep","edep*1000.0",proton+direct);
    hits->Project("ex","edep*1000.0",xtalk);

    hits->Project("rconv","sqrt(shpos.y^2+shpos.x^2)",conv+direct);
//    hits->Project("rdio","sqrt(shpos.y^2+shpos.x^2)",dio+direct);
//    hits->Project("rneut","sqrt(shpos.y^2+shpos.x^2)",neutron);
//    hits->Project("rphot","sqrt(shpos.y^2+shpos.x^2)",photon);
    hits->Project("rdelta","sqrt(shpos.y^2+shpos.x^2)",bkge+direct);
    hits->Project("rp","sqrt(shpos.y^2+shpos.x^2)",proton+direct);
//    hits->Project("rx","sqrt(shpos.y^2+shpos.x^2)",xtalk);
    /*    
	  hits->Project("nconv","n200",conv);
	  hits->Project("ndio","n200",bkge);
	  hits->Project("ndelta","n200",bkgo);
	  hits->Project("np","n200",proton);

	  hits->Project("zconv","shpos.z",conv);
	  hits->Project("zdio","shpos.z",bkge);
	  hits->Project("zdelta","shpos.z",bkgo);
	  hits->Project("zp","shpos.z",proton);
     */
    TCanvas* bcan = new TCanvas("bcan","background",1000,800);
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
//    edelta->GetXaxis()->SetRangeUser(-0.005,0.04);
    edelta->SetMinimum(1);
    double maxv = edelta->GetMaximum();
    maxv = max(maxv, ep->GetMaximum());
    maxv = max(maxv, ex->GetMaximum());
    edelta->SetMaximum(2*maxv);
    edelta->Draw();
    ep->Draw("same");
    econv->Draw("same");
//    edio->Draw("same");
    ex->Draw("same");
//    eneut->Draw("same");
//    ephot->Draw("same");

    TBox* eselbox = new TBox(0.0,edelta->GetMinimum(),3.5,edelta->GetMaximum());
//    Int_t tYellow = TColor::GetColorTransparent(kYellow,0.3);
    eselbox->SetFillColor(kYellow);
    eselbox->SetFillStyle(3004);
//    eselbox->SetLineStyle(3);
    eselbox->Draw();
//    TLine* ecut_t = new TLine(0.003,0.0,0.003,edelta->GetMaximum());
//    ecut_t->SetLineColor(kBlack);
//    ecut_t->SetLineStyle(2);
//    ecut_t->SetLineWidth(2);
//    TLine* ecut_l = new TLine(0.0,0.0,0.0,edelta->GetMaximum());
//    ecut_l->SetLineColor(kBlack);
//    ecut_l->SetLineStyle(2);
//    ecut_l->SetLineWidth(2);
//    ecut_t->Draw();
//    ecut_l->Draw();

    TLegend* leg2 = new TLegend(0.55,0.7,0.9,0.9);
    leg2->AddEntry(rconv,"CE Induced","l");
    leg2->AddEntry(rdelta,"Background e Induced","l");
//    leg2->AddEntry(rdio,"DIO Electrons","l");
//    leg2->AddEntry(rneut,"Neutrons","l");
//    leg2->AddEntry(rphot,"Photons","l");
    leg2->AddEntry(rp,"Proton Induced","l");
    leg2->AddEntry(ex,"X-talk","l");
    leg2->AddEntry(eselbox,"Track Reconstruction Selection","F");
//    leg3->AddEntry(ecut_l,"Loose cut","l");
    leg2->Draw();

    int istart = econv->FindFixBin(0.0);
    int istop = econv->FindFixBin(0.003);
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
//    rdio->Draw("same");
//    rx->Draw("same");
//    rneut->Draw("same");
 //   rphot->Draw("same");

    TBox* rselbox = new TBox(395.0,rdelta->GetMinimum(),650,rdelta->GetMaximum());
    rselbox->SetFillColor(kYellow);
    rselbox->SetFillStyle(3004);
//    rselbox->SetLineStyle(3);
    rselbox->Draw();

//    TLine* rmin_t = new TLine(395,0.0,395,rdelta->GetMaximum());
//    rmin_t->SetLineColor(kBlack);
//    rmin_t->SetLineStyle(2);
//    rmin_t->SetLineWidth(2);
//    TLine* rmin_l = new TLine(390,0.0,390,rdelta->GetMaximum());
//    rmin_l->SetLineColor(kBlack);
//    rmin_l->SetLineStyle(3);
//    rmin_l->SetLineWidth(2);

//    TLine* rmax_t = new TLine(650,0.0,650,rdelta->GetMaximum());
//    rmax_t->SetLineColor(kBlack);
//    rmax_t->SetLineStyle(2);
//    rmax_t->SetLineWidth(2);
//    TLine* rmax_l = new TLine(650,0.0,650,rdelta->GetMaximum());
//    rmax_l->SetLineColor(kBlack);
//    rmax_l->SetLineStyle(3);
//    rmax_l->SetLineWidth(2);
//    rmin_t->Draw();
//    rmin_l->Draw();
//    rmax_t->Draw();
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
  } else if(spage=="origin"){

    TCut timecut("time>500");

    THStack* origin = new THStack("origin","Reco Hit Time by Generator Particle;Hit Time (ns);Hits/event/ns");
    TH1F* dtime = new TH1F("dtime","DIO Reco Hit Time",150,250,1750);
    TH1F* gtime = new TH1F("gtime","Photon Reco Hit Time",150,250,1750);
    TH1F* pptime = new TH1F("pptime","Primary Proton Reco Hit Time",150,250,1750);
    TH1F* stptime = new TH1F("stptime","Stopping Target Proton Reco Hit Time",150,250,1750);
    TH1F* ntime = new TH1F("ntime","Neutron Reco Hit Time",150,250,1750);
    TH1F* mtime = new TH1F("mtime","OOT Muon Reco Hit Time",150,250,1750);
    TH1F* cetime = new TH1F("cetime","CE Reco Hit Time",150,250,1750);
    dtime->SetFillColor(kOrange);
    gtime->SetFillColor(kBlack);
    pptime->SetFillColor(kBlue);
    stptime->SetFillColor(kGreen);
    ntime->SetFillColor(kCyan);
    mtime->SetFillColor(kYellow);
    cetime->SetFillColor(kRed);

    double scale = 0.1/nevents;
    hits->Project("cetime","time",conv+timecut);
    cetime->Scale(scale);
    origin->Add(cetime);
    hits->Project("dtime","time",dioorigin+timecut);
    dtime->Scale(scale);
    origin->Add(dtime);
    hits->Project("mtime","time",ootmuonorigin+timecut);
    mtime->Scale(scale);
    origin->Add(mtime);
    hits->Project("gtime","time",porigin+timecut);
    gtime->Scale(scale);
    origin->Add(gtime);
    hits->Project("stptime","time",stpprotonorigin+timecut);
    stptime->Scale(scale);
    origin->Add(stptime);
    hits->Project("pptime","time",pprotonorigin+timecut);
    pptime->Scale(scale);
    origin->Add(pptime);
    hits->Project("ntime","time",norigin+timecut);
    ntime->Scale(scale);
    origin->Add(ntime);
 
    double total = dtime->Integral() + pptime->Integral() + stptime->Integral() + gtime->Integral() + ntime->Integral() + mtime->Integral() + cetime->Integral() ;
 
    cout << "DIO integral = " << dtime->Integral()
      << " Primary Proton inegral = " << pptime->Integral()
      << " ST Proton inegral = " << stptime->Integral()
      << " Photon inegral = " << gtime->Integral()
      << " Neutron inegral = " << ntime->Integral()
      << " OOT muon inegral = " << mtime->Integral()
      << " CE inegral = " << cetime->Integral()
      << " total = " << total << endl;

    TCanvas* ocan = new TCanvas("ocan","ocan",800,800);
    origin->Draw("h");
    TLegend* tleg = new TLegend(.6,.6,.9,.9);
    char title[50];
    snprintf(title,50,"Neutron, #int=%4.0f",ntime->Integral()*10.0);
    tleg->AddEntry(ntime,title,"F");
    snprintf(title,50,"Primary Proton, #int=%4.0f",pptime->Integral()*10.0);
    tleg->AddEntry(pptime,title,"F");
    snprintf(title,50,"Stopping Target Proton, #int=%4.0f",stptime->Integral()*10.0);
    tleg->AddEntry(stptime,title,"F");
    snprintf(title,50,"Photon, #int=%4.0f",gtime->Integral()*10.0);
    tleg->AddEntry(gtime,title,"F");
    snprintf(title,50,"OOT Muon, #int=%4.0f",mtime->Integral()*10.0);
    tleg->AddEntry(mtime,title,"F");
    snprintf(title,50,"DIO, #int=%4.0f",dtime->Integral()*10.0);
    tleg->AddEntry(dtime,title,"F");
    snprintf(title,50,"CE, #int=%4.0f",cetime->Integral()*10.0);
    tleg->AddEntry(cetime,title,"F");
    snprintf(title,50,"Total #int=%4.0f",total*10.0);
    tleg->AddEntry(dtime,title,"");
    tleg->Draw();
  } else if(spage=="sorigin"){

    THStack* sorigin = new THStack("sorigin","Selected Hit Time by Generator Particle;Hit Time (ns);Hits/event/ns");
    TH1F* dtime = new TH1F("dtime","DIO Reco Hit Time",150,250,1750);
    TH1F* gtime = new TH1F("gtime","Photon Reco Hit Time",150,250,1750);
    TH1F* pptime = new TH1F("pptime","Primary Proton Reco Hit Time",150,250,1750);
    TH1F* stptime = new TH1F("stptime","Stopping Target Proton Reco Hit Time",150,250,1750);
    TH1F* ntime = new TH1F("ntime","Neutron Reco Hit Time",150,250,1750);
    TH1F* mtime = new TH1F("mtime","OOT Muon Reco Hit Time",150,250,1750);
    TH1F* cetime = new TH1F("cetime","CE Reco Hit Time",150,250,1750);
    dtime->SetFillColor(kOrange);
    gtime->SetFillColor(kBlack);
    pptime->SetFillColor(kBlue);
    stptime->SetFillColor(kGreen);
    ntime->SetFillColor(kCyan);
    mtime->SetFillColor(kYellow);
    cetime->SetFillColor(kRed);

    double scale = 0.1/nevents;
    hits->Project("cetime","mctime",conv+hitsel);
    cetime->Scale(scale);
    sorigin->Add(cetime);
     hits->Project("dtime","mctime",dioorigin+hitsel);
    dtime->Scale(scale);
    sorigin->Add(dtime);
    hits->Project("mtime","mctime",ootmuonorigin+hitsel);
    mtime->Scale(scale);
    sorigin->Add(mtime);
    hits->Project("gtime","mctime",porigin+hitsel);
    gtime->Scale(scale);
    sorigin->Add(gtime);
    hits->Project("stptime","mctime",stpprotonorigin+hitsel);
    stptime->Scale(scale);
    sorigin->Add(stptime);
    hits->Project("pptime","mctime",pprotonorigin+hitsel);
    pptime->Scale(scale);
    sorigin->Add(pptime);
    hits->Project("ntime","mctime",norigin+hitsel);
    ntime->Scale(scale);
    sorigin->Add(ntime);


    double total = dtime->Integral() + pptime->Integral() + stptime->Integral() + gtime->Integral() + ntime->Integral() + mtime->Integral() + cetime->Integral();

    cout << "DIO integral = " << dtime->Integral()
      << " Primary Proton inegral = " << pptime->Integral()
      << " ST Proton inegral = " << stptime->Integral()
      << " Photon inegral = " << gtime->Integral()
      << " Neutron inegral = " << ntime->Integral()
      << " OOT muon inegral = " << mtime->Integral()
      << " CE inegral = " << cetime->Integral()
      << " total = " << total << endl;

    TCanvas* socan = new TCanvas("socan","socan",800,800);
    sorigin->Draw("h");
    TLegend* tleg = new TLegend(.6,.6,.9,.9);
    char title[50];
    snprintf(title,50,"Neutron, #int=%4.0f",ntime->Integral()*10.0);
    tleg->AddEntry(ntime,title,"F");
    snprintf(title,50,"Primary Proton, #int=%4.0f",pptime->Integral()*10.0);
    tleg->AddEntry(pptime,title,"F");
    snprintf(title,50,"Stopping Target Proton, #int=%4.0f",stptime->Integral()*10.0);
    tleg->AddEntry(stptime,title,"F");
    snprintf(title,50,"Photon, #int=%4.0f",gtime->Integral()*10.0);
    tleg->AddEntry(gtime,title,"F");
    snprintf(title,50,"OOT Muon, #int=%4.0f",mtime->Integral()*10.0);
    tleg->AddEntry(mtime,title,"F");
    snprintf(title,50,"DIO, #int=%4.0f",dtime->Integral()*10.0);
    tleg->AddEntry(dtime,title,"F");
    snprintf(title,50,"CE, #int=%4.0f",cetime->Integral()*10.0);
    tleg->AddEntry(cetime,title,"F");
    snprintf(title,50,"Total #int=%4.0f",total*10.0);
    tleg->AddEntry(dtime,title,"");
    tleg->Draw();
  } else if(spage=="hitsel"){
    TH2F* hsel = new TH2F("hsel","Hit Selection",5,-0.5,4.5,7,-0.5,6.5);
    TAxis* yax = hsel->GetYaxis();
    unsigned ibin(1);
    yax->SetBinLabel(ibin++,"All");
    yax->SetBinLabel(ibin++,"Radius");
    yax->SetBinLabel(ibin++,"Time");
    yax->SetBinLabel(ibin++,"Energy");
    yax->SetBinLabel(ibin++,"Isolated");
    yax->SetBinLabel(ibin++,"Cluster");
    yax->SetBinLabel(ibin++,"Good");
    TAxis* xax = hsel->GetXaxis();
    ibin = 1;
    xax->SetBinLabel(ibin++,"Proton");
    xax->SetBinLabel(ibin++,"LowEE");
    xax->SetBinLabel(ibin++,"DIO");
    xax->SetBinLabel(ibin++,"Other");
    xax->SetBinLabel(ibin++,"CE");
    // first, get normalization
    TH1F* myhp = new TH1F("myhp","My hit particle",5,-0.5,4.5);
    TH1F* myhpg = new TH1F("myhpg","My hit particle",5,-0.5,4.5);
    xax = myhp->GetXaxis();
    ibin = 1;
    xax->SetBinLabel(ibin++,"Proton");
    xax->SetBinLabel(ibin++,"LowEE");
    xax->SetBinLabel(ibin++,"DIO");
    xax->SetBinLabel(ibin++,"Other");
    xax->SetBinLabel(ibin++,"CE");
    hits->Project("myhp","myhpart(mcpdg,mcgen,mcproc)");
    hits->Project("myhpg","myhpart(mcpdg,mcgen,mcproc)",hitsel);
    myhpg->SetFillColor(kGreen);
    // now loop over selections
    std::vector<TCut> selcuts = {"","rsel","tsel","esel","isolated","delta",hitsel};
    for(size_t icut=0;icut< selcuts.size();++icut){
      char val[100];
      cout << "Projecting cut " << selcuts[icut] << endl;
      snprintf(val,100,"%lu:myhpart(mcpdg,mcgen,mcproc)",icut);
      hits->Project("+hsel",val,selcuts[icut]);
    }
// normalize by row
    for(int ibin=1;ibin <= myhp->GetXaxis()->GetNbins();++ibin){
      double norm = 1.0/myhp->GetBinContent(ibin);
      cout << "Normalization for " << hsel->GetYaxis()->GetBinLabel(ibin)  << " = " << norm << endl;
      for(int jbin=1;jbin <= hsel->GetYaxis()->GetNbins(); ++jbin) {
	double val =hsel->GetBinContent(ibin,jbin);
	cout << "value for ibin " << ibin <<" jbin " << jbin << " val " << val << endl;
	hsel->SetBinContent(ibin,jbin,val*norm);
      }
    }
    TCanvas* hscan = new TCanvas("hscan","hscan",800,800);
    hscan->Divide(1,2);
    hscan->cd(1);
    hsel->Draw("box");
    hscan->cd(2);
    myhp->Draw();
    myhpg->Draw("same");
  
  } 
}
