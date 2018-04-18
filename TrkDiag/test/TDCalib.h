//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Dec 27 08:50:36 2017 by ROOT version 6.10/06
// from TTree shdiag/strawhit diagnostics
// found on file: /data/SHDM.root
//////////////////////////////////////////////////////////

#ifndef TDCalib_h
#define TDCalib_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TCut.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include <vector>

// Header file for the classes stored in the TTree if any.

class TDCalib {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
   static constexpr Int_t kMaxshpos = 1;
   static constexpr Int_t kMaxmcshpos = 1;
   static constexpr Int_t kMaxmcopos = 1;
   static constexpr Int_t kMaxmcpopos = 1;
   static constexpr Int_t kMaxmcgpos = 1;

   // Declaration of leaf types
   Int_t           eventid;
 //CLHEP::Hep3Vector *shpos_;
   Double_t        shpos_dx;
   Double_t        shpos_dy;
   Double_t        shpos_dz;
   Float_t         shlen;
   Float_t         slen;
   Float_t         edep;
   Float_t         time_tcal;
   Float_t         time_thv;
   Float_t         tot_totcal;
   Float_t         tot_tothv;
   Float_t         rho;
   Int_t           plane;
   Int_t           panel;
   Int_t           layer;
   Int_t           straw;
   Int_t           esel;
   Int_t           rsel;
   Int_t           tsel;
   Int_t           bkgclust;
   Int_t           bkg;
   Int_t           stereo;
   Int_t           tdiv;
   Int_t           strawxtalk;
   Int_t           elecxtalk;
   Int_t           isolated;
   Int_t           calosel;
   Float_t         pdist;
   Float_t         pperp;
   Float_t         pmom;
   Float_t         wres;
   Float_t         tres;
   Float_t         shchisq;
   Float_t         shdt;
   Float_t         shdist;
 //CLHEP::Hep3Vector *mcshpos_;
   Double_t        mcshpos_dx;
   Double_t        mcshpos_dy;
   Double_t        mcshpos_dz;
 //CLHEP::Hep3Vector *mcopos_;
   Double_t        mcopos_dx;
   Double_t        mcopos_dy;
   Double_t        mcopos_dz;
 //CLHEP::Hep3Vector *mcpopos_;
   Double_t        mcpopos_dx;
   Double_t        mcpopos_dy;
   Double_t        mcpopos_dz;
   Float_t         mcct_mcctcal;
   Float_t         mcct_mccthv;
   Float_t         mcoe;
   Float_t         mcom;
   Float_t         mcpoe;
   Float_t         mcpom;
   Float_t         mcshlen;
   Float_t         mcshd;
   Float_t         mcplen;
   Float_t         mcedep;
   Float_t         mcetrig;
   Int_t           mcnsteps;
   Int_t           mcpdg;
   Int_t           mcgen;
   Int_t           mcproc;
   Float_t         mcsptime;
   Float_t         mcwt_mcwtcal;
   Float_t         mcwt_mcwthv;
   Int_t           mcppdg;
   Int_t           mcpproc;
   Double_t        mcptime;
   Int_t           mcgid;
   Int_t           mcgpdg;
   Float_t         mcge;
   Float_t         mcgt;
 //CLHEP::Hep3Vector *mcgpos_;
   Double_t        mcgpos_dx;
   Double_t        mcgpos_dy;
   Double_t        mcgpos_dz;
   Char_t          mcxtalk;

   // List of branches
   TBranch        *b_eventid;   //!
   TBranch        *b_shpos_dx;   //!
   TBranch        *b_shpos_dy;   //!
   TBranch        *b_shpos_dz;   //!
   TBranch        *b_shlen;   //!
   TBranch        *b_slen;   //!
   TBranch        *b_edep;   //!
   TBranch        *b_time;   //!
   TBranch        *b_tot;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_plane;   //!
   TBranch        *b_panel;   //!
   TBranch        *b_layer;   //!
   TBranch        *b_straw;   //!
   TBranch        *b_esel;   //!
   TBranch        *b_rsel;   //!
   TBranch        *b_tsel;   //!
   TBranch        *b_bkgclust;   //!
   TBranch        *b_bkg;   //!
   TBranch        *b_stereo;   //!
   TBranch        *b_tdiv;   //!
   TBranch        *b_strawxtalk;   //!
   TBranch        *b_elecxtalk;   //!
   TBranch        *b_isolated;   //!
   TBranch        *b_calosel;   //!
   TBranch        *b_pdist;   //!
   TBranch        *b_pperp;   //!
   TBranch        *b_pmom;   //!
   TBranch        *b_wres;   //!
   TBranch        *b_tres;   //!
   TBranch        *b_shchisq;   //!
   TBranch        *b_shdt;   //!
   TBranch        *b_shdist;   //!
   TBranch        *b_mcshpos_dx;   //!
   TBranch        *b_mcshpos_dy;   //!
   TBranch        *b_mcshpos_dz;   //!
   TBranch        *b_mcopos_dx;   //!
   TBranch        *b_mcopos_dy;   //!
   TBranch        *b_mcopos_dz;   //!
   TBranch        *b_mcpopos_dx;   //!
   TBranch        *b_mcpopos_dy;   //!
   TBranch        *b_mcpopos_dz;   //!
   TBranch        *b_mcct;   //!
   TBranch        *b_mcoe;   //!
   TBranch        *b_mcom;   //!
   TBranch        *b_mcpoe;   //!
   TBranch        *b_mcpom;   //!
   TBranch        *b_mcshlen;   //!
   TBranch        *b_mcshd;   //!
   TBranch        *b_mcplen;   //!
   TBranch        *b_mcedep;   //!
   TBranch        *b_mcetrig;   //!
   TBranch        *b_mcnsteps;   //!
   TBranch        *b_mcpdg;   //!
   TBranch        *b_mcgen;   //!
   TBranch        *b_mcproc;   //!
   TBranch        *b_mcsptime;   //!
   TBranch        *b_mcwt;   //!
   TBranch        *b_mcppdg;   //!
   TBranch        *b_mcpproc;   //!
   TBranch        *b_mcptime;   //!
   TBranch        *b_mcgid;   //!
   TBranch        *b_mcgpdg;   //!
   TBranch        *b_mcge;   //!
   TBranch        *b_mcgt;   //!
   TBranch        *b_mcgpos_dx;   //!
   TBranch        *b_mcgpos_dy;   //!
   TBranch        *b_mcgpos_dz;   //!
   TBranch        *b_mcxtalk;   //!

   TDCalib(float emin=0.2,float emax=7.0, float ebin=0.2, TTree* ttree=0);
   virtual ~TDCalib();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

  float _emin, _emax, _ebin;
  float _central;

  std::vector<TCanvas*> _cans;
  std::vector<TH1F*> _evec;
  std::vector<TH2F*> _tdevec;
  std::vector<TH2F*> _resvec;
  TH1F *_ceedep, *_pedep, *_oedep;
  TH2F *_pulledep, *_pullwlen, *_pulltedep, *_pulltwlen, *_pullsfrac, *_pulltsfrac;
  TH1F *_cepull, *_ppull, *_ceres, *_pres;
  std::vector<TProfile*> _tdpevec;
  std::vector<float> _emean, _eerr;
  std::vector<float> _smean, _serr;
  std::vector<float> _cres, _creserr;
  std::vector<float> _sres, _sreserr;
  TGraphErrors *_cresg, *_sresg;
 
  void MakeHists();
  void SaveCans(const char* suffix);
  void SaveData(const char* datafile);
};

#endif

#ifdef TDCalib_cxx
TDCalib::TDCalib( float emin, float emax, float ebin, TTree *tree) : fChain(0),
  _emin(emin), _emax(emax) ,_ebin(ebin), _central(65.0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/data/SHDM.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/data/SHDM.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("/data/SHDM.root:/SHD");
      dir->GetObject("shdiag",tree);

   }
   Init(tree);
}

TDCalib::~TDCalib()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t TDCalib::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t TDCalib::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void TDCalib::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("eventid", &eventid, &b_eventid);
   fChain->SetBranchAddress("shpos.dx", &shpos_dx, &b_shpos_dx);
   fChain->SetBranchAddress("shpos.dy", &shpos_dy, &b_shpos_dy);
   fChain->SetBranchAddress("shpos.dz", &shpos_dz, &b_shpos_dz);
   fChain->SetBranchAddress("shlen", &shlen, &b_shlen);
   fChain->SetBranchAddress("slen", &slen, &b_slen);
   fChain->SetBranchAddress("edep", &edep, &b_edep);
   fChain->SetBranchAddress("time", &time_tcal, &b_time);
   fChain->SetBranchAddress("tot", &tot_totcal, &b_tot);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("plane", &plane, &b_plane);
   fChain->SetBranchAddress("panel", &panel, &b_panel);
   fChain->SetBranchAddress("layer", &layer, &b_layer);
   fChain->SetBranchAddress("straw", &straw, &b_straw);
   fChain->SetBranchAddress("esel", &esel, &b_esel);
   fChain->SetBranchAddress("rsel", &rsel, &b_rsel);
   fChain->SetBranchAddress("tsel", &tsel, &b_tsel);
   fChain->SetBranchAddress("bkgclust", &bkgclust, &b_bkgclust);
   fChain->SetBranchAddress("bkg", &bkg, &b_bkg);
   fChain->SetBranchAddress("stereo", &stereo, &b_stereo);
   fChain->SetBranchAddress("tdiv", &tdiv, &b_tdiv);
   fChain->SetBranchAddress("strawxtalk", &strawxtalk, &b_strawxtalk);
   fChain->SetBranchAddress("elecxtalk", &elecxtalk, &b_elecxtalk);
   fChain->SetBranchAddress("isolated", &isolated, &b_isolated);
   fChain->SetBranchAddress("calosel", &calosel, &b_calosel);
   fChain->SetBranchAddress("pdist", &pdist, &b_pdist);
   fChain->SetBranchAddress("pperp", &pperp, &b_pperp);
   fChain->SetBranchAddress("pmom", &pmom, &b_pmom);
   fChain->SetBranchAddress("wres", &wres, &b_wres);
   fChain->SetBranchAddress("tres", &tres, &b_tres);
   fChain->SetBranchAddress("shchisq", &shchisq, &b_shchisq);
   fChain->SetBranchAddress("shdt", &shdt, &b_shdt);
   fChain->SetBranchAddress("shdist", &shdist, &b_shdist);
   fChain->SetBranchAddress("mcshpos.dx", &mcshpos_dx, &b_mcshpos_dx);
   fChain->SetBranchAddress("mcshpos.dy", &mcshpos_dy, &b_mcshpos_dy);
   fChain->SetBranchAddress("mcshpos.dz", &mcshpos_dz, &b_mcshpos_dz);
   fChain->SetBranchAddress("mcopos.dx", &mcopos_dx, &b_mcopos_dx);
   fChain->SetBranchAddress("mcopos.dy", &mcopos_dy, &b_mcopos_dy);
   fChain->SetBranchAddress("mcopos.dz", &mcopos_dz, &b_mcopos_dz);
   fChain->SetBranchAddress("mcpopos.dx", &mcpopos_dx, &b_mcpopos_dx);
   fChain->SetBranchAddress("mcpopos.dy", &mcpopos_dy, &b_mcpopos_dy);
   fChain->SetBranchAddress("mcpopos.dz", &mcpopos_dz, &b_mcpopos_dz);
   fChain->SetBranchAddress("mcct", &mcct_mcctcal, &b_mcct);
   fChain->SetBranchAddress("mcoe", &mcoe, &b_mcoe);
   fChain->SetBranchAddress("mcom", &mcom, &b_mcom);
   fChain->SetBranchAddress("mcpoe", &mcpoe, &b_mcpoe);
   fChain->SetBranchAddress("mcpom", &mcpom, &b_mcpom);
   fChain->SetBranchAddress("mcshlen", &mcshlen, &b_mcshlen);
   fChain->SetBranchAddress("mcshd", &mcshd, &b_mcshd);
   fChain->SetBranchAddress("mcplen", &mcplen, &b_mcplen);
   fChain->SetBranchAddress("mcedep", &mcedep, &b_mcedep);
   fChain->SetBranchAddress("mcetrig", &mcetrig, &b_mcetrig);
   fChain->SetBranchAddress("mcnsteps", &mcnsteps, &b_mcnsteps);
   fChain->SetBranchAddress("mcpdg", &mcpdg, &b_mcpdg);
   fChain->SetBranchAddress("mcgen", &mcgen, &b_mcgen);
   fChain->SetBranchAddress("mcproc", &mcproc, &b_mcproc);
   fChain->SetBranchAddress("mcsptime", &mcsptime, &b_mcsptime);
   fChain->SetBranchAddress("mcwt", &mcwt_mcwtcal, &b_mcwt);
   fChain->SetBranchAddress("mcppdg", &mcppdg, &b_mcppdg);
   fChain->SetBranchAddress("mcpproc", &mcpproc, &b_mcpproc);
   fChain->SetBranchAddress("mcptime", &mcptime, &b_mcptime);
   fChain->SetBranchAddress("mcgid", &mcgid, &b_mcgid);
   fChain->SetBranchAddress("mcgpdg", &mcgpdg, &b_mcgpdg);
   fChain->SetBranchAddress("mcge", &mcge, &b_mcge);
   fChain->SetBranchAddress("mcgt", &mcgt, &b_mcgt);
   fChain->SetBranchAddress("mcgpos.dx", &mcgpos_dx, &b_mcgpos_dx);
   fChain->SetBranchAddress("mcgpos.dy", &mcgpos_dy, &b_mcgpos_dy);
   fChain->SetBranchAddress("mcgpos.dz", &mcgpos_dz, &b_mcgpos_dz);
   fChain->SetBranchAddress("mcxtalk", &mcxtalk, &b_mcxtalk);
   Notify();

  MakeHists();

}

Bool_t TDCalib::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void TDCalib::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
#endif // #ifdef TDCalib_cxx
