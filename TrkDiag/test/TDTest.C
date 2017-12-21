#include "TTree.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TDirectory.h"
#include "TStyle.h"
#include "TCut.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <string>
#include <vector>
#include <iostream>
using std::cout;
using std::endl;

Double_t ALine(Double_t *x, Double_t *par) {
  if(x[0] > par[0])
    return par[1];
  else
    return par[1]+(x[0]-par[0])*par[2];
}
 
void TDTest(TTree* sh,const char* page="edep") {
  std::string spage(page);
  TCut conv("mcproc==56&&mcgen==2");
  TCut proton("mcpdg==2212");
  if(spage == "edep"){
    TH1F* cedep = new TH1F("cedep","Ce Straw Energy Deposit;Edep (KeV)",50,0,8.0);
    TH1F* pedep = new TH1F("pedep","Proton Straw Energy Deposit;Edep (KeV)",50,0,8.0);
    cedep->SetStats(0);
    pedep->SetStats(0);
    cedep->SetLineColor(kRed);
    pedep->SetLineColor(kBlue);
    sh->Project("cedep","edep*1000.0",conv);
    sh->Project("pedep","edep*1000.0",proton);
    cedep->Scale(pedep->GetEntries()/cedep->GetEntries());
    TCanvas* edep = new TCanvas("edep","edep",400,400);
    cedep->Draw("Hist");
    pedep->Draw("same");
    TLegend* leg = new TLegend(0.6,0.7,0.9,0.9);
    leg->AddEntry(cedep,"Conversion electron","l");
    leg->AddEntry(pedep,"Proton","l");
    leg->Draw();

  } else if(spage == "edlen"){
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(0);
    TH2F* edslen = new TH2F("edslen","Ce Reco Hit Energy vs Straw Half Length;Straw Half Length (mm);EDep (KeV)",40,200.0,600.0,50,0,6.0);
    TProfile* edslenp = new TProfile("edslenp","Ce Reco Hit Energy vs Straw Half Length;Straw Half Length (mm);EDep (KeV)",40,200.0,600.0,0.0,3.0,"S");
    edslen->SetStats(0);
    sh->Project("edslen","1000.0*edep:slen",conv);
    sh->Project("edslenp","1000.0*edep:slen",conv);

    edslen->FitSlicesY(0,0,-1,20);

    TCanvas* ed = new TCanvas("ed","ed",600,600);
    TH1D *edslen_1 = (TH1D*)gDirectory->Get("edslen_1");
    edslen_1->SetLineColor(kRed);
    edslen_1->Fit("pol1");
    edslen_1->SetStats(1);
    edslen->Draw("colorz");
    edslen_1->Draw("sames");
    edslenp->Draw("same");
  } else if(spage == "TDlen"){
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(0);
    TH2F* ctdlen = new TH2F("ctdlen","Ce TDRes vs reco wire position;reco wire position (mm);reco - mc wire position (mm)",20,0.0,400.0,50,-100.0,100.0);
    TH2F* ptdlen = new TH2F("ptdlen","Proton TDRes vs reco wire position;reco wire position (mm);reco - mc wire position (mm)",20,0.0,400.0,50,-80.0,80.0);
    sh->Project("ctdlen","shlen-mcshlen:abs(shlen)",conv);
    sh->Project("ptdlen","shlen-mcshlen:abs(shlen)",proton);
    ctdlen->FitSlicesY(0,0,-1);
    ptdlen->FitSlicesY(0,0,-1);
    TCanvas* tdlcan = new TCanvas("tdlcan","tdlcan",700,700);
    tdlcan->Divide(2,2);
    tdlcan->cd(1);
    ctdlen->Draw("colorz");
    tdlcan->cd(2);
    ptdlen->Draw("colorz");
    tdlcan->cd(3);
    TH1D *ctdlen_2 = (TH1D*)gDirectory->Get("ctdlen_2");
    ctdlen_2->SetLineColor(kRed);
    ctdlen_2->Fit("pol3");
    ctdlen_2->SetStats(1);
    tdlcan->cd(4);
    TH1D *ptdlen_2 = (TH1D*)gDirectory->Get("ptdlen_2");
    ptdlen_2->SetLineColor(kRed);
    ptdlen_2->Fit("pol3");
    ptdlen_2->SetStats(1);
    TH2F* ctdlenp = new TH2F("ctdlenp","Ce TD pull vs reco wire position;reco wire position (mm);TD pull",15,0.0,450.0,50,-6.0,6.0);
    TH2F* ptdlenp = new TH2F("ptdlenp","Proton TD pull vs reco wire position;reco wire position (mm);TD pull",15,0.0,450.0,50,-6.0,6.0);
    sh->Project("ctdlenp","(shlen-mcshlen)/wres:abs(shlen)",conv);
    sh->Project("ptdlenp","(shlen-mcshlen)/wres:abs(shlen)",proton);
    ctdlenp->FitSlicesY(0,0,-1);
    ptdlenp->FitSlicesY(0,0,-1);
    TCanvas* tdlpcan = new TCanvas("tdlpcan","tdlpcan",700,700);
    tdlpcan->Divide(2,2);
    tdlpcan->cd(1);
    ctdlenp->Draw("colorz");
    tdlpcan->cd(2);
    ptdlenp->Draw("colorz");
    tdlpcan->cd(3);
    TH1D *ctdlenp_2 = (TH1D*)gDirectory->Get("ctdlenp_2");
    ctdlenp_2->SetLineColor(kRed);
    ctdlenp_2->Fit("pol3");
    ctdlenp_2->SetStats(1);
    tdlpcan->cd(4);
    TH1D *ptdlenp_2 = (TH1D*)gDirectory->Get("ptdlenp_2");
    ptdlenp_2->SetLineColor(kRed);
    ptdlenp_2->Fit("pol3");
    ptdlenp_2->SetStats(1);
  } else if(spage == "TDREdep"){
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(0);
    std::vector<TH2F*> tdrvec;
    char cut[50];
    char name[50];
    char title[80];
    unsigned nplots(0);
    float emin(0.0);
    float ebin(0.5);
    float emax = emin+ebin;
    while(emax < 7.0){
      ++nplots;
      snprintf(cut,50,"1000.0*edep>%3.1f&&1000.0*edep<%3.1f",emin,emax);
      TCut slcut(cut);
      snprintf(name,50,"tdr%3.1f",emin);
      snprintf(title,80,"TD pull vs reco wire position, %3.1fKeV < edep < %3.1fKeV;reco wire position (mm);TD pull",emin,emax);
      TH2F* hist = new TH2F(name,title,10,0.0,450.0,50,-6.0,6.0);
      tdrvec.push_back(hist);
      sh->Project(name,"(shlen-mcshlen)/wres:abs(shlen)",cut);
      emin += ebin;
      emax += ebin;
    }
    TCanvas* tdrcan = new TCanvas("tdrcan","tdrcan",700,700);
    unsigned nx = (unsigned)ceil(sqrt(nplots));
    unsigned ny = ceil(nplots/nx);
    tdrcan->Divide(nx,ny);
    for(unsigned iplot=0;iplot<nplots;++iplot){
      tdrcan->cd(iplot+1);
      tdrvec[iplot]->Draw("colorz");
    }
    TCanvas* tdrfcan = new TCanvas("tdrfcan","tdrfcan",700,700);
    tdrfcan->Divide(nx,ny);
    for(unsigned iplot=0;iplot<nplots;++iplot){
      tdrfcan->cd(iplot+1);
      tdrvec[iplot]->FitSlicesY(0,0,-1);
      std::string fname(tdrvec[iplot]->GetName());
      fname += "_2";
      TH1D *tdr_2 = (TH1D*)gDirectory->Get(fname.c_str());
      tdr_2->SetLineColor(kRed);
      tdr_2->Fit("pol3");
      tdr_2->SetStats(1);
    }
  } else if(spage == "TDEdep"){
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(0);
    std::vector<TH1F*> evec;
    std::vector<TH1F*> rvec;
    std::vector<TH2F*> tdevec;
    std::vector<TProfile*> tdepvec;
    char cut[50];
    char name[50];
    char pname[50];
    char ename[50];
    char rname[50];
    char title[80];
    char rtitle[80];
    char resname[80];
    unsigned nplots(0);
    float emin(0.2);
    float ebin(0.2);
//    float ebin(1.0);
    float emax = emin+ebin;
    while(emax < 7.0){
      ++nplots;
      snprintf(cut,50,"1000.0*edep>%3.1f&&1000.0*edep<%3.1f&&abs(shlen)<450.0",emin,emax);
      TCut slcut(cut);
      snprintf(name,50,"td%3.1f",emin);
      snprintf(title,80,"True wire position vs #Delta t, %3.1fKeV < edep < %3.1fKeV",emin,emax);
      TH2F* hist = new TH2F(name,title,50,-10.0,10.0,40,-600.0,600.0);
      tdevec.push_back(hist);
      sh->Project(name,"mcshlen:tcal-thv",cut);
      snprintf(pname,50,"tdp%3.1f",emin);
      TProfile* prof = new TProfile(pname,title,50,-10.0,10.0,-600.0,600.0);
      tdepvec.push_back(prof);
      sh->Project(pname,"mcshlen:tcal-thv",cut);
      snprintf(ename,50,"edep%3.1f",emin);
      TH1F* edep = new TH1F(ename,title,100,0.0,7.0);
      evec.push_back(edep);
      sh->Project(ename,"1000.0*edep",cut);
      emin += ebin;
      emax += ebin;
    }
    TCanvas* tdecan = new TCanvas("tdecan","tdecan",700,700);
    unsigned nx = (unsigned)ceil(sqrt(nplots));
    unsigned ny = ceil(nplots/nx);
    TF1* line = new TF1("line","[0]+[1]*x",-10.0,10.0);
    line->SetParameter(0,0.0);
    line->SetParameter(1,100.0);
    tdecan->Divide(nx,ny);
    std::vector<float> emean, eerr, rmean, rerr;
    std::vector<float> smean, serr;
    emin = 0.0;
    emax = emin+ebin;
    for(unsigned iplot=0;iplot<nplots;++iplot){
      tdecan->cd(iplot+1);
      tdevec[iplot]->Draw("colorz");
      if(tdepvec[iplot]->GetEntries() > 20){
	TFitResultPtr fitr = tdepvec[iplot]->Fit(line,"QS","sames",-3.5,3.5);
	double slope = fitr->Parameter(1);
	//      tdpvec[iplot]->Draw("sames");
	snprintf(rname,50,"res%3.1f",emin);
	snprintf(rtitle,80,"Wire Distance Resolution, %3.1fKeV < edep < %3.1fKeV",emin,emax);
	snprintf(resname,80,"(tcal-thv)*%8.3g-mcshlen",slope);
	TH1F* res = new TH1F(rname,rtitle,150,-300.0,300.0);
	rvec.push_back(res);
	sh->Project(rname,resname,cut);
	smean.push_back(slope);
	serr.push_back(fitr->ParError(1));
	emean.push_back(evec[iplot]->GetMean());
	eerr.push_back(evec[iplot]->GetRMS());
      }
      emin += ebin;
      emax += ebin;
    }

    TCanvas* edcan = new TCanvas("edcan","edcan",700,700);
    edcan->Divide(nx,ny);
    for(unsigned iplot=0;iplot<nplots;++iplot){
      edcan->cd(iplot+1);
      evec[iplot]->Draw();
    }
  
    TCanvas* tdrcan = new TCanvas("tdrcan","tdrcan",700,700);
    tdrcan->Divide(nx,ny);
    TF1* ngaus = new TF1("ngaus","0.398942280401*[0]*exp( 0.5*( -((x-[1])/[2])^2))/[2]",-10.0,10.0);
    for(unsigned iplot=0;iplot<nplots;++iplot){
      tdrcan->cd(iplot+1);
      ngaus->SetParameters(rvec[iplot]->GetEntries(),0.0,rvec[iplot]->GetRMS());
      ngaus->SetRange(-1.5*rvec[iplot]->GetRMS(),1.5*rvec[iplot]->GetRMS());
      TFitResultPtr gfit = rvec[iplot]->Fit(ngaus,"QS","",-2.0*rvec[iplot]->GetRMS(),2.0*rvec[iplot]->GetRMS());
      rmean.push_back(gfit->Parameter(2));
      rerr.push_back(gfit->ParError(2));
    }
    // fit slope vs edep
    TGraphErrors* slopes = new TGraphErrors(nplots,emean.data(),smean.data(),eerr.data(),serr.data());
    slopes->SetTitle("TD Slope vs EDep;EDep (KeV);Slope (mm/ns)");
    TF1* sfit = new TF1("ALine",ALine,0.0,10.0,3);
    sfit->SetParameters(3.0,130.0,21.0);
    slopes->Fit(sfit,"Q");
    // fit resolution vs edep
    TGraphErrors* tdres = new TGraphErrors(nplots,emean.data(),rmean.data(),eerr.data(),rerr.data());
    tdres->SetTitle("TD Resolution vs EDep;EDep (KeV);Resolution #sigma (mm)");
    TF1* rfit = new TF1("ALine",ALine,0.0,10.0,3);
    rfit->SetParameters(3.0,20.0,-21.0);
//    rfit->FixParameter(0,sfit->GetParameter(0));
    tdres->Fit(rfit,"Q");

    TCanvas* fits = new TCanvas("fits","fits",700,700);
    fits->Divide(1,2);
    fits->cd(1);
    slopes->Draw("ALP");
    fits->cd(2);
    tdres->Draw("ALP");

// pull plots
    char comps[120];
    snprintf(comps,120,"min(%5.2f,%5.2f+(1000.0*edep-%5.2f)*%5.2f)*(tcal-thv):mcshlen",
    sfit->GetParameter(1), sfit->GetParameter(1),
    sfit->GetParameter(0), sfit->GetParameter(2)); 
    cout << "comp string = " << comps << endl;
    TH2F* ncomp = new TH2F("ncomp","new Ce TD Reco wire distance vs true;True TD distance (mm);Reco wire distance (mm)",50,-600,600,50,-600,600);
    TH2F* ocomp = new TH2F("ocomp","old Ce TD Reco wire distance vs true;True TD distance (mm);Reco wire distance (mm)",50,-600,600,50,-600,600);
    sh->Project("ncomp",comps,conv);
    sh->Project("ocomp","shlen:mcshlen",conv);

    char ress[120];
    snprintf(ress,120,"min(%5.2f,%5.2f+(1000.0*edep-%5.2f)*%5.2f)*(tcal-thv)-mcshlen",
    sfit->GetParameter(1), sfit->GetParameter(1),
    sfit->GetParameter(0), sfit->GetParameter(2)); 
    TH1F* nres = new TH1F("nres","new Ce TD Wire Distance Resolution; Resolution (mm)",50,-300,300);
    TH1F* ores = new TH1F("ores","old Ce TD Wire Distance Resolution; Resolution (mm)",50,-300,300);
    sh->Project("nres",ress,conv);
    sh->Project("ores","shlen-mcshlen",conv);
    TF1* diag = new TF1("diag","x",-600,600);
  
    char pulls[120];
    snprintf(pulls,120,"(min(%5.2f,%5.2f+(1000.0*edep-%5.2f)*%5.2f)*(tcal-thv)-mcshlen)/max(%5.2f,%5.2f+(1000.0*edep-%5.2f)*%5.2f)",
    sfit->GetParameter(1), sfit->GetParameter(1),
    sfit->GetParameter(0), sfit->GetParameter(2),
    rfit->GetParameter(1), rfit->GetParameter(1),
    rfit->GetParameter(0), rfit->GetParameter(2)); 

    TH1F* npull = new TH1F("npull","new TD Ce Pull;Pull",100,-10,10);
    TH1F* opull = new TH1F("opull","old TD Ce Pull;Pull",100,-10,10);
    sh->Project("npull",pulls,conv);
    sh->Project("opull","(shlen-mcshlen)/wres",conv);
    gStyle->SetOptStat(1111);
  
    TCanvas* occan = new TCanvas("occan","occan",700,700);
    occan->Divide(2,2);
    occan->cd(1);
    ocomp->Draw("colorz");
    diag->Draw("same");
    occan->cd(2);
    ores->Fit("gaus","","",-2.0*ores->GetRMS(),2.0*ores->GetRMS());
    occan->cd(3);
    opull->Fit("gaus","","",-2.0*opull->GetRMS(),2.0*opull->GetRMS());

    TCanvas* nccan = new TCanvas("nccan","nccan",700,700);
    nccan->Divide(2,2);
    nccan->cd(1);
    ncomp->Draw("colorz");
    diag->Draw("same");
    nccan->cd(2);
    nres->Fit("gaus","","",-2.0*nres->GetRMS(),2.0*nres->GetRMS());
    nccan->cd(3);
    npull->Fit("gaus","","",-2.0*npull->GetRMS(),2.0*npull->GetRMS());
  } else if(spage == "reco"){
    TH2F* tdce = new TH2F("tdce","Ce Reco TD wire distance vs true;True TD distance (mm);Reco wire distance (mm)",50,-650,650,50,-650,650);
    TH2F* tdp = new TH2F("tdp","Proton TD Reco wire distance vs true;True TD distance (mm);Reco wire distance (mm)",50,-650,650,50,-650,650);
    tdce->SetStats(0);
    tdp->SetStats(0);
    sh->Project("tdce","shlen:mcshlen",conv);
    sh->Project("tdp","shlen:mcshlen",proton);
    TCanvas* rcan = new TCanvas("rcan","rcan",800,400);
    rcan->Divide(2,1);
    rcan->cd(1);
    tdce->Draw("colorz");
    rcan->cd(2);
    tdp->Draw("colorz");
  }
}
