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
  enum mp{ beamproton=0, proton=1, photon=2, neutron=3, electron=4, ootmuon=5, deuteron=6, other=7};
  if(part==2212&&gen==16)return beamproton;
  if(part==2212&&gen==28)return proton;
  if(part==1000010020&&gen==28)return deuteron;
  if(part==22) return photon;
  if(part==11)return electron;
  if(part==2112)return neutron;
  if(part==13)return ootmuon;
  return other;
}


Int_t mypart(Int_t part) {
  enum mp{ proton=0, electron=1, photon=2, neutron=3, positron=4, muon=5, deuteron=6, other=7};
  if(part==2212)return proton;
  if(part==1000010020)return deuteron;
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
  enum mp{ primary=0, compt=1, bkg=2, gamconv=3, photo=4, had=5, xt=6, other=7};
  if(xtalk)return xt;
  if(proc==56)return primary;
  if(proc==12)return compt;
  if(proc==17 || proc==21) return bkg;
  if(proc==13)return gamconv;
  if(proc==40)return photo;
  if(proc==58)return had;
  return other;
}

Int_t myhpart(Int_t mcpdg,Int_t mcgen, Int_t mcproc){
  enum hpart{Hadron=0,LowEe,Mudecay,CE,Other};
  if(mcgen==2)return CE;
  if(mcpdg==11&&(mcgen==28||mcproc==14||mcproc==114))return Mudecay;
  if(abs(mcpdg)==11)return LowEe;
  if(mcpdg>2000)return Hadron;
  return Other;
}

void StrawHitTest (TTree* hits, const char* page="bcan",unsigned nevents=1000 ) {

  TString spage(page);
//  gStyle->SetOptStat(0);
//  TCut conv("mcpdg==11&&mcgen==2&&mcmom>100.0");
  TCut conv("mcpdg==11&&mcgen==2");
  TCut oele("abs(mcpdg)==11&&mcgen!=2");
//  TCut dio("mcpdg==11&&mcgen==6"); MC truth bug
  TCut dio("mcpdg==11&&(mcproc==14||mcproc==56&&mcoe<90)");
  TCut bkg("mcpdg==11&&mcgen<0&&mcproc==17");
  TCut pconv("mcpdg==11&&mcgen<0&&mcproc==11");
  TCut compt("mcpdg==11&&mcgen<0&&mcproc==12");
  TCut proton("mcpdg==2212");
  TCut deuteron("mcpdg==1000010020");
  TCut hadron("mcpdg>2000");
  TCut photon("mcpdg==22");
  TCut neutron("mcpdg==2112");
  TCut bkge("abs(mcpdg)==11&&mcppdg==22");
  TCut bkgo("abs(mcpdg)!=11&&mcpdg!=2212");
  TCut xtalk("mcxtalk");
  TCut direct("!mcxtalk");
  TCut opart=!(conv||oele||proton);

  TCut pproton("mcppdg==2212");
  TCut pdeuteron("mcppdg==1000010020");
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
  TCut stpdeuteronorigin("mcgpdg==1000010020&&mcgid==28");
  TCut pprotonorigin("mcgpdg==2212&&mcgid==16");
  TCut flashstraw("plane>33&&straw>=90");

  TCut hitsel("esel&&rsel&&tsel&&(!bkg)");

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

    TH1F* nuconv = new TH1F("nuconv","N steps in StrawHit",6,-0.5,5.5);
    TH1F* nudio = new TH1F("nudio","N steps in StrawHit",6,-0.5,5.5);
    TH1F* nudel = new TH1F("nudel","N steps in StrawHit",6,-0.5,5.5);
    TH1F* nup = new TH1F("nup","N steps in StrawHit",6,-0.5,5.5);
    nuconv->SetLineColor(kRed);
    nudio->SetLineColor(kGreen);
    nudel->SetLineColor(kCyan);
    nup->SetLineColor(kBlue);
    TH2F* eve = new TH2F("eve","StrawHit Edep vs #Sigma MC Edep;MeV;MeV",100,0,0.05,100,0,0.05);
    eve->SetStats(0);

    hits->Project("pdgid","mcpdg");
    hits->Project("gide","mcgen","mcpdg==11");
    hits->Project("gidp","mcgen",proton);
    hits->Project("pide","mcproc",bkg);
    hits->Project("pidp","mcproc",proton);
    hits->Project("nuconv","mcnsteps",conv);
    hits->Project("nudio","mcnsteps",dio);
    hits->Project("nudel","mcnsteps",bkg);
    hits->Project("nup","mcnsteps",proton);

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

  } else if(spage=="particle"){
//    gStyle->SetOptStat(111111);
    THStack* estack = new THStack("edep","Reco Hit Energy by Particle;Deposited Energy (KeV)");
    TH1F* econv = new TH1F("econv","Straw Hit Energy;Deposited Energy (KeV)",200,-1.0,12.0);
    estack->Add(econv);
    TH1F* emu = new TH1F("emu","Straw Hit Energy;Deposited Energy (KeV)",200,-1.0,12.0);
    estack->Add(emu);
    TH1F* egam = new TH1F("egam","Straw Hit Energy;Deposited Energy (KeV)",200,-1.0,12.0);
    estack->Add(egam);
    TH1F* ehad = new TH1F("ehad","Straw Hit Energy;Deposited Energy (KeV)",200,-1.0,12.0);
    estack->Add(ehad);
    econv->SetFillColor(kRed);
    emu->SetFillColor(kCyan);
    egam->SetFillColor(kGreen);
    ehad->SetFillColor(kMagenta);
//    econv->SetStats(0);
    emu->SetStats(0);
    egam->SetStats(0);
    ehad->SetStats(0);

    THStack* rstack = new THStack("rho","Reco Hit Radius by Particle;Transverse Radius (mm)");
    TH1F* rconv = new TH1F("rconv","Straw Hit Radius;Transverse Radius (mm)",100,360,720);
    rstack->Add(rconv);
    TH1F* rmu = new TH1F("rmu","Straw Hit Radius;Transverse Radius (mm)",100,360,720);
    rstack->Add(rmu);
    TH1F* rgam = new TH1F("rgam","Straw Hit Radius;Transverse Radius (mm)",100,360,720);
    rstack->Add(rgam);
    TH1F* rhad = new TH1F("rhad","Straw Hit Radius;Transverse Radius (mm)",100,360,720);
    rstack->Add(rhad);
//    TH1F* rx = new TH1F("rx","Straw Hit Radius;Transverse Radius (mm)",100,360,720);
    rconv->SetFillColor(kRed);
    rmu->SetFillColor(kCyan);
    rgam->SetFillColor(kGreen);
    rhad->SetFillColor(kMagenta);
    hits->Project("econv","edep*1000.0",conv+direct);
    hits->Project("emu","edep*1000.0",dio+direct);
    hits->Project("egam","edep*1000.0",bkge+direct);
    hits->Project("ehad","edep*1000.0",hadron+direct);

    hits->Project("rconv","sqrt(shpos.dy^2+shpos.dx^2)",conv+direct);
    hits->Project("rmu","sqrt(shpos.dy^2+shpos.dx^2)",dio+direct);
    hits->Project("rgam","sqrt(shpos.dy^2+shpos.dx^2)",bkge+direct);
    hits->Project("rhad","sqrt(shpos.dy^2+shpos.dx^2)",hadron+direct);
    
    TCanvas* bcan = new TCanvas("bcan","background",1000,800);
    bcan->Divide(1,2);
    bcan->cd(1);
//    gPad->SetLogy();
    estack->Draw();
    econv->Draw("same");
    TLegend* leg2 = new TLegend(0.55,0.7,0.9,0.9);
    leg2->AddEntry(econv,"Conversion electron","f");
    leg2->AddEntry(emu,"#mu#rightarrowe#nu#nu","f");
    leg2->AddEntry(egam,"Low-E e","f");
    leg2->AddEntry(ehad,"Hadron","f");
    leg2->Draw();

    bcan->cd(2);
 //   gPad->SetLogy();
    rstack->Draw();
  } else if (spage == "ccan") {

    TCanvas* ccan = new TCanvas("ccan","cleaned hits",1200,800);
    TCut clean = goodpeak+TCut("tight>0&&bkg==0");
        
	  TH1F* gid = new TH1F("gid","Generator code",21,-1.5,19.5);
	  TH1F* gidc = new TH1F("gidc","Generator code",21,-1.5,19.5);

	  TH1F* rres = new TH1F("rres","StrawHit Radius resolution;mm",100,-200,200);
   
    TH1F* pres = new TH1F("pres","StrawHit #phi resolution;rad",100,-0.5,0.5);
        

	  gid->SetLineColor(kBlue);
	  gidc->SetLineColor(kRed);
	  hits->Project("gid","mcgen",clean);
	  hits->Project("gidc","mcgen",conv);

	  hits->Project("rres","sqrt(shpos.dy^2+shpos.dx^2)-sqrt(mcshpos.dy^2+mcshpos.dx^2)");
    
    hits->Project("pres","atan2(shpos.dy,shpos.dx)-atan2(mcshpos.dy,mcshpos.dx)");

    /*
       TLegend* leg3 = new TLegend(0.4,0.75,0.7,0.9);
       leg3->AddEntry(gidc,"Conv. Electron hits","l");
       leg3->AddEntry(gid,"Selected hits","l");
       leg3->Draw();
     */
//    gStyle->SetOptStat(0);
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
    TH1F* tbkg = new TH1F("tbkg","time WRT peak;nsec",101,-80,80);
    TH1F* tp = new TH1F("tp","time WRT peak;nsec",101,-80,80);
    tconv->SetLineColor(kRed);
    tdio->SetLineColor(kGreen);
    tbkg->SetLineColor(kCyan);
    tp->SetLineColor(kBlue);

    hits->Project("tconv","time-tpeak",conv+goodevt+clean);
    hits->Project("tdio","time-tpeak",bkge+goodevt+clean);
    hits->Project("tbkg","time-tpeak",bkgo+goodevt+clean);
    hits->Project("tp","time-tpeak",proton+goodevt+clean);

    TLegend* leg4 = new TLegend(0.65,0.6,0.9,0.95);
    leg4->AddEntry(tconv,"Conv. Electrons","l");
    leg4->AddEntry(tdio,"Bkg Electrons","l");
    leg4->AddEntry(tbkg,"Photons","l");
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
    tbkg->Draw("same");
    tp->Draw("same");
    leg4->Draw();
    tmin->Draw();
    tmax->Draw();

    // compute efficiency and purity
    Float_t ncon = gidc->GetEntries();
    Float_t nsel = gid->GetEntries();
    Float_t nsbkg = gid->GetBinContent(1);
    Float_t nscon = gid->GetBinContent(4);
    Float_t nsdio = gid->GetBinContent(8);
    Float_t nsprot = gid->GetBinContent(14);
    cout << "selection efficiency = " << nscon/ncon << endl;
    cout << "selection purity = " << nscon/nsel << endl;
    cout << "bkg fraction bkgs = " << nsbkg/nsel << " dio " << nsdio/nsel << " proton " << nsprot/nsel << endl;
  } else if(spage=="origin"){

    TCut timecut("time>450");

    THStack* origin = new THStack("origin","Reco Hit Time by Generator Particle;Hit Time (ns);Hits/event/ns");
    TH1F* dtime = new TH1F("dtime","DIO Reco Hit Time",300,250,1750);
    TH1F* gtime = new TH1F("gtime","Photon Reco Hit Time",300,250,1750);
    TH1F* pptime = new TH1F("pptime","Primary Proton Reco Hit Time",300,250,1750);
    TH1F* stptime = new TH1F("stptime","Muon Daughter Proton Reco Hit Time",300,250,1750);
    TH1F* stdtime = new TH1F("stdtime","Muon Daughter Detueron Reco Hit Time",300,250,1750);
    TH1F* ntime = new TH1F("ntime","Neutron Reco Hit Time",300,250,1750);
    TH1F* mtime = new TH1F("mtime","OOT Muon Reco Hit Time",300,250,1750);
    TH1F* cetime = new TH1F("cetime","CE Reco Hit Time",300,250,1750);
    dtime->SetFillColor(kOrange);
    gtime->SetFillColor(kBlack);
    pptime->SetFillColor(kBlue);
    stptime->SetFillColor(kGreen);
    stdtime->SetFillColor(kMagenta);
    ntime->SetFillColor(kCyan);
    mtime->SetFillColor(kYellow);
    cetime->SetFillColor(kRed);

    double scale = 1.0/(dtime->GetBinWidth(1)*nevents);
    hits->Project("cetime","time",conv+timecut);
    cetime->Scale(scale);
    origin->Add(cetime);
    hits->Project("dtime","time",dioorigin+timecut);
    dtime->Scale(scale);
    origin->Add(dtime);
    hits->Project("stdtime","time",stpdeuteronorigin+timecut);
    stdtime->Scale(scale);
    origin->Add(stdtime);
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
      << " ST Deuteron inegral = " << stdtime->Integral()
      << " CE inegral = " << cetime->Integral()
      << " total = " << total << endl;

    TCanvas* ocan = new TCanvas("ocan","ocan",800,800);
    origin->Draw("h");
    TLegend* tleg = new TLegend(.6,.6,.9,.9);
    char title[50];
    double factor = dtime->GetBinWidth(1);
    snprintf(title,50,"Neutron, #int=%4.0f",ntime->Integral()*factor);
    tleg->AddEntry(ntime,title,"F");
    snprintf(title,50,"Primary Proton, #int=%4.0f",pptime->Integral()*factor);
    tleg->AddEntry(pptime,title,"F");
    snprintf(title,50,"Stopping Target Proton, #int=%4.0f",stptime->Integral()*factor);
    tleg->AddEntry(stptime,title,"F");
    snprintf(title,50,"Photon, #int=%4.0f",gtime->Integral()*factor);
    tleg->AddEntry(gtime,title,"F");
    snprintf(title,50,"OOT Muon, #int=%4.0f",mtime->Integral()*factor);
    tleg->AddEntry(mtime,title,"F");
    snprintf(title,50,"Deuteron, #int=%4.0f",stdtime->Integral()*factor);
    tleg->AddEntry(stdtime,title,"F");
    snprintf(title,50,"DIO, #int=%4.0f",dtime->Integral()*factor);
    tleg->AddEntry(dtime,title,"F");
    snprintf(title,50,"CE, #int=%4.0f",cetime->Integral()*factor);
    tleg->AddEntry(cetime,title,"F");
    snprintf(title,50,"Total #int=%4.0f",total*factor);
    tleg->AddEntry(dtime,title,"");
    tleg->Draw();
  } else if(spage=="sorigin"){

    THStack* sorigin = new THStack("sorigin","Selected Hit Time by Generator Particle;Hit Time (ns);Hits/event/ns");
    TH1F* dtime = new TH1F("dtime","DIO Reco Hit Time",150,250,1750);
    TH1F* gtime = new TH1F("gtime","Photon Reco Hit Time",150,250,1750);
    TH1F* pptime = new TH1F("pptime","Primary Proton Reco Hit Time",150,250,1750);
    TH1F* stptime = new TH1F("stptime","Stopping Target Proton Reco Hit Time",150,250,1750);
    TH1F* stdtime = new TH1F("stdtime","Stopping Target Deuteron Reco Hit Time",150,250,1750);
    TH1F* ntime = new TH1F("ntime","Neutron Reco Hit Time",150,250,1750);
    TH1F* mtime = new TH1F("mtime","OOT Muon Reco Hit Time",150,250,1750);
    TH1F* cetime = new TH1F("cetime","CE Reco Hit Time",150,250,1750);
    dtime->SetFillColor(kOrange);
    gtime->SetFillColor(kBlack);
    pptime->SetFillColor(kBlue);
    stptime->SetFillColor(kGreen);
    stdtime->SetFillColor(kMagenta);
    ntime->SetFillColor(kCyan);
    mtime->SetFillColor(kYellow);
    cetime->SetFillColor(kRed);

    double scale = 1.0/(dtime->GetBinWidth(1)*nevents);
    hits->Project("cetime","time",conv+hitsel);
    cetime->Scale(scale);
    sorigin->Add(cetime);
     hits->Project("dtime","time",dioorigin+hitsel);
    dtime->Scale(scale);
    sorigin->Add(dtime);
    hits->Project("stdtime","time",stpdeuteronorigin+hitsel);
    stdtime->Scale(scale);
    sorigin->Add(stdtime);
    hits->Project("mtime","time",ootmuonorigin+hitsel);
    mtime->Scale(scale);
    sorigin->Add(mtime);
    hits->Project("gtime","time",porigin+hitsel);
    gtime->Scale(scale);
    sorigin->Add(gtime);
    hits->Project("stptime","time",stpprotonorigin+hitsel);
    stptime->Scale(scale);
    sorigin->Add(stptime);
    hits->Project("pptime","time",pprotonorigin+hitsel);
    pptime->Scale(scale);
    sorigin->Add(pptime);
    hits->Project("ntime","time",norigin+hitsel);
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
    double factor = dtime->GetBinWidth(1);
    snprintf(title,50,"Neutron, #int=%4.0f",ntime->Integral()*factor);
    tleg->AddEntry(ntime,title,"F");
    snprintf(title,50,"Primary Proton, #int=%4.0f",pptime->Integral()*factor);
    tleg->AddEntry(pptime,title,"F");
    snprintf(title,50,"Stopping Target Proton, #int=%4.0f",stptime->Integral()*factor);
    tleg->AddEntry(stptime,title,"F");
    snprintf(title,50,"Photon, #int=%4.0f",gtime->Integral()*factor);
    tleg->AddEntry(gtime,title,"F");
    snprintf(title,50,"OOT Muon, #int=%4.0f",mtime->Integral()*factor);
    tleg->AddEntry(mtime,title,"F");
    snprintf(title,50,"Stopping Target Deuteron, #int=%4.0f",stdtime->Integral()*factor);
    tleg->AddEntry(stdtime,title,"F");
    snprintf(title,50,"DIO, #int=%4.0f",dtime->Integral()*factor);
    tleg->AddEntry(dtime,title,"F");
    snprintf(title,50,"CE, #int=%4.0f",cetime->Integral()*factor);
    tleg->AddEntry(cetime,title,"F");
    snprintf(title,50,"Total #int=%4.0f",total*factor);
    tleg->AddEntry(dtime,title,"");
    tleg->Draw();
  } else if(spage=="hitsel"){
    TH2F* hsel = new TH2F("hsel","Hit Selection;Producing Particle;Cut efficiency (%)",5,-0.5,4.5,4,-0.5,3.5);
    TAxis* yax = hsel->GetYaxis();
    unsigned ibin(1);
    yax->SetBinLabel(ibin++,"Hit Time");
    yax->SetBinLabel(ibin++,"Hit Energy");
    yax->SetBinLabel(ibin++,"Hit Radius");
    yax->SetBinLabel(ibin++,"Bkg Hit");
    TAxis* xax = hsel->GetXaxis();
    ibin = 1;
    xax->SetBinLabel(ibin++,"Hadron");
    xax->SetBinLabel(ibin++,"LowEe");
    xax->SetBinLabel(ibin++,"#mu#rightarrowe#nu#nu");
    xax->SetBinLabel(ibin++,"Ce");
    xax->SetBinLabel(ibin++,"Other");
    // first, get normalization
    TH1F* myhp = new TH1F("myhp","Hit Rate;Producing Particle;Hits/event",5,-0.5,4.5);
    TH1F* myhpg = new TH1F("myhpg","Hit Rate;Producing Particle;Hits/event",5,-0.5,4.5);
    xax = myhp->GetXaxis();
    ibin = 1;
    xax->SetBinLabel(ibin++,"Hadron");
    xax->SetBinLabel(ibin++,"LowEe");
    xax->SetBinLabel(ibin++,"#mu#rightarrowe#nu#nu");
    xax->SetBinLabel(ibin++,"Ce");
    xax->SetBinLabel(ibin++,"Other");
    hits->Project("myhp","myhpart(mcpdg,mcgen,mcproc)","tsel");
    hits->Project("myhpg","myhpart(mcpdg,mcgen,mcproc)",hitsel);
    double pscale(1.0/nevents);
    myhp->Scale(pscale);
    myhpg->Scale(pscale);
    myhpg->SetFillColor(kGreen);
    // now loop over selections
    std::vector<TCut> selcuts = {"tsel","tsel&&esel","tsel&&esel&&rsel",hitsel};
    for(size_t icut=0;icut< selcuts.size();++icut){
      char val[100];
      cout << "Projecting cut " << selcuts[icut] << endl;
      snprintf(val,100,"%lu:myhpart(mcpdg,mcgen,mcproc)",icut);
      hits->Project("+hsel",val,selcuts[icut]);
    }
// normalize by row
    for(int ibin=1;ibin <= myhp->GetXaxis()->GetNbins();++ibin){
      double norm = 100.0/hsel->GetBinContent(ibin,1);
      cout << "Normalization for " << hsel->GetYaxis()->GetBinLabel(ibin)  << " = " << norm << endl;
      for(int jbin=1;jbin <= hsel->GetYaxis()->GetNbins(); ++jbin) {
	double val =hsel->GetBinContent(ibin,jbin);
	cout << "value for ibin " << ibin <<" jbin " << jbin << " val " << val << endl;
	hsel->SetBinContent(ibin,jbin,val*norm);
      }
    }
    TCanvas* hscan = new TCanvas("hscan","hscan",750,750);
    hscan->Divide(1,2);
    hscan->cd(1);
    hsel->Draw("boxtext0");
    hscan->cd(2);
    TLegend* leg = new TLegend(0.6,0.7,0.8,0.9);
    leg->AddEntry(myhp,"All Hits","l");
    leg->AddEntry(myhpg,"Selected Hits","f");
    myhp->Draw("histtext0");
    myhpg->Draw("histtext90same");
    leg->Draw();
  
  } else if(spage == "tot") {

    TH2F* ptot = new TH2F("ptot","Proton TOT vs MC Transverse Drift Distance;True Drift Distance (mm);TOT (ns)",50,0,2.5,16,0,64);
    TH2F* etot = new TH2F("etot","Electron TOT vs MC Transverse Drift Distance;True Drift Distance (mm);TOT (ns)",50,0,2.5,16,0,64);
    ptot->SetStats(0);
    etot->SetStats(0);
    hits->Project("ptot","0.5*(totcal+tothv):abs(mcshd)","mcpdg==2212");
    hits->Project("etot","0.5*(totcal+tothv):abs(mcshd)","mcpdg==11&&mcproc==56&&mcoe>100");
    TCanvas* totcan = new TCanvas("totcan","TOT can",800,600);
    totcan->Divide(2,1);
    totcan->cd(1);
    etot->Draw("colorz");
    totcan->cd(2);
    ptot->Draw("colorz");

  } else if(spage == "td") {
    TH1F* tdres = new TH1F("tdres","Time Division Resolution",100,-400,400);
    hits->Project("tdres","shlen-mcshlen","mcgen==2");
    TCanvas* tdcan = new TCanvas("tdcan","tdcan",600,600);
    tdres->Fit("gaus");
  } else if(spage == "flash") {
    TH1F* tfhit = new TH1F("tfhit","Flash Hit Time",500,0.0,250.0);
    TH1F* tfmc = new TH1F("tfmc","Flash Hit Time",500,0.0,250.0);
    hits->Project("tfhit","tcal",pprotonorigin+flashstraw);
    hits->Project("tfmc","mcsptime",pprotonorigin+flashstraw);
    tfhit->SetLineColor(kRed);
    tfmc->SetLineColor(kBlue);
    TCanvas* fcan = new TCanvas("fcan","fcan",600,600);
    fcan->cd(0);
    tfhit->Draw();
    tfmc->Draw("same");
  }
}
