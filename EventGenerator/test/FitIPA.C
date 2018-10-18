#include "TFile.h"
#include "TStyle.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TBranch.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TCut.h"
#include "Rtypes.h"
struct XYZT {
  Float_t x, y, z, time;
};

Double_t EMG(Double_t *x, Double_t *par) {
  static const double invsq2 = 1.0/sqrt(2.0);
  double ls2 = par[3]*par[2]*par[2];
  Double_t result = 0.5*par[0]*par[3]*TMath::Exp(0.5*par[3]*(2.0*(par[1]-x[0])+ls2));
  result *= TMath::Erfc(invsq2*(par[1]-x[0] + ls2 )/par[2]);
  return result;
}

class FitIPA {
  public:
  FitIPA(const char* filename="/cvmfs/mu2e.opensciencegrid.org/DataFiles/mergedMuonStops/nts.mu2e.DS-OOTstops.MDC2018a.001002_00000000.root",
      const char* dirname="outOfTargetDumper");
  void Fit();
  void Fill();
  void CreateIPAStops (const char* filename="IPAstops.root",unsigned nrows=10000);
  TCanvas* ipacan, *tcan;
  TFile *ipafile;
  TDirectory * ootd;
  TFile* newfile;
  TDirectory* ipadir;
  TTree* stops, *newstops;
  TF1* sine;
  TF1* zline;
  TF1* tgaus;
  TH1F* radius;
  TH1F* zpos;
  TH1F* phi;
  TH1F* ipamutime;
  TH2F* xy;
  TCut loosez;
  TCut ipaz;
  TCut ipar;
  double xoff;
};


FitIPA::FitIPA(const char* filename,const char* dirname) :
  loosez("z>6900&&z<8000"),
  ipaz("z>6920&&z<7920"),
  xoff(3904)
{
  gStyle->SetOptFit(111);
  gStyle->SetOptStat(110011);
  ipafile = TFile::Open(filename);
  ootd = (TDirectory*)ipafile->Get(dirname);
  stops = (TTree*)ootd->Get("stops");
  sine = new TF1("sine","[0]*sin(x+[1])+[2]",-3.14159,3.14159);
  sine->SetParameters(4.0,3.0,5.0);
  zline = new TF1("zline","[0]+[1]*x",6920,7920);
  zline->SetParameters(-35,0.005);
  tgaus = new TF1("tgaus",EMG,100,400,4);
  tgaus->SetParameters(500,146,16.5,0.036);
  xy = new TH2F("xy","OOT mu stop;x (mm);y (mm)",50,-350,350,50,-350,350);
  xy->SetStats(0);
  radius = new TH1F("radius","IPA mu stop radius;r (mm)",30,299.0,300.5);
  zpos = new TH1F("zpos","IPA mu stop z;z (mm)",25,6920,7920);
  phi = new TH1F("phi","IPA mu stop azimuth;#phi (rad)",25,-3.14159,3.14159);
  char cutstring[100];
  snprintf(cutstring,100,"sqrt(y^2+(x+%f)^2)-299.75 < 0.25",xoff);
  ipar = TCut (cutstring);
}

void FitIPA::Fill() {
  char pstring[100];
  snprintf(pstring,100,"sqrt(y^2+(x+%f)^2)",xoff);
  stops->Project("radius",pstring,ipaz);
  stops->Project("zpos","z",ipar&&loosez);
  snprintf(pstring,100,"atan2(y,x+%f)",xoff);
  stops->Project("phi",pstring,ipar&&ipaz);
}

void FitIPA::Fit() {
  ipacan = new TCanvas("ipacan","ipacan",600,600);
  ipacan->Divide(2,2);
  ipacan->cd(1);
  xy->Draw();
  char pstring[100];
  snprintf(pstring,100,"y:x+%f",xoff);
  stops->Draw(pstring,loosez,"same");
  ipacan->cd(2);
  radius->Draw();
  ipacan->cd(3);
//  zpos->SetMaximum(25.0);
  zpos->Fit(zline);
  ipacan->cd(4);
  sine->SetParameter(0,phi->GetEntries()/30.0);
  phi->Fit(sine);
  ipamutime = new TH1F("ipamutime","IPA mu stop time",40,0,400);
  stops->Project("ipamutime","time",ipaz&&ipar);
  tcan = new TCanvas("tcan","tcan",400,400);
  tgaus->SetParameter(0,10*ipamutime->GetEntries());
  ipamutime->Fit(tgaus);
}
void FitIPA::CreateIPAStops (const char* filename,unsigned nrows) {
  newfile = new TFile(filename,"RECREATE");
  ipadir = newfile->mkdir("IPAstops","IPAstops");
  newfile->cd("IPAstops");
  newstops = new TTree("stops","Stopped particles ntuple");
  XYZT xyzt;
  newstops->Branch("stops",&xyzt,"x/F:y/F:z/F:time/F");
  TRandom3 myrand;
  for(unsigned irow=0;irow < nrows;++irow){
    double phi = sine->GetRandom();
    double radius = myrand.Uniform(299.5,300.0);
    double z = zline->GetRandom();
    double t = tgaus->GetRandom();
    xyzt.x = radius*cos(phi) - xoff;
    xyzt.y = radius*sin(phi);
    xyzt.z = z;
    xyzt.time = t;
    newstops->Fill();
  }
  cout << "TTree has " <<  newstops->GetEntries() << " entries" << endl;
  newstops->Write();
  newfile->Write();
  newfile->Close();
}
