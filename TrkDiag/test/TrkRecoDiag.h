//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Dec 15 20:34:36 2016 by ROOT version 6.08/02
// from TChain TrkRecoDiag/trdiag/
//////////////////////////////////////////////////////////

#ifndef TrkRecoDiag_h
#define TrkRecoDiag_h

#include <TROOT.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <vector>

// Header file for the classes stored in the TTree if any.

class TrkRecoDiag {
public :
   TTree          *_chain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
   const Int_t kMaxhsf = 1;
   const Int_t kMaxhsh = 1;
   const Int_t kMaxksf = 1;
   const Int_t kMaxksh = 1;
   const Int_t kMaxkff = 1;
   const Int_t kMaxkfh = 1;
   const Int_t kMaxmch = 1;

   // Declaration of leaf types
   Int_t           iev;
   Float_t         tct0;
   Float_t         tct0err;
   Int_t           tcn;
 //mu2e::BitMap<mu2e::TrkFitFlagDetail> *hsf_;
   UInt_t          hsf__value;
 //mu2e::RobustHelix *hsh_;
   Float_t         hsh__rcent;
   Float_t         hsh__fcent;
   Float_t         hsh__radius;
   Float_t         hsh__lambda;
   Float_t         hsh__fz0;
   Float_t         hst0;
   Float_t         hst0err;
   Int_t           hsn;
   Int_t           hsna;
 //mu2e::BitMap<mu2e::TrkFitFlagDetail> *ksf_;
   UInt_t          ksf__value;
 //mu2e::HelixVal  *ksh_;
   Float_t         ksh__pars[5];
   Float_t         kst0;
   Float_t         kst0err;
   Float_t         ksm;
   Float_t         ksmerr;
   Int_t           ksn;
   Int_t           ksna;
 //mu2e::BitMap<mu2e::TrkFitFlagDetail> *kff_;
   UInt_t          kff__value;
 //mu2e::HelixVal  *kfh_;
   Float_t         kfh__pars[5];
   Float_t         kft0;
   Float_t         kft0err;
   Float_t         kfm;
   Float_t         kfmerr;
   Int_t           kfn;
   Int_t           kfna;
   Float_t         beamwt;
   UInt_t          ndigitot;
   Float_t         mcgenmom;
   Float_t         mcentmom;
   Float_t         mcentpz;
   Float_t         mcentt0;
   Float_t         mcmidmom;
   Float_t         mcmidpz;
   Float_t         mcmidt0;
   Float_t         mcxitmom;
   Float_t         mcxitpz;
   Float_t         mcxitt0;
 //mu2e::RobustHelix *mch_;
   Float_t         mch__rcent;
   Float_t         mch__fcent;
   Float_t         mch__radius;
   Float_t         mch__lambda;
   Float_t         mch__fz0;
   Int_t           pdg;
   Int_t           gen;
   Int_t           proc;
   UInt_t          ndigi;
   UInt_t          nkfprimary;
   UInt_t          nksprimary;
   UInt_t          nhsprimary;
   UInt_t          ntcprimary;
   Int_t           kfnp;
   Int_t           kfnap;
   Int_t           ksnp;
   Int_t           ksnap;
   Int_t           hsnp;
   Int_t           hsnap;
   Int_t           tcnp;

   // List of branches
   TBranch        *b_iev;   //!
   TBranch        *b_tct0;   //!
   TBranch        *b_tct0err;   //!
   TBranch        *b_tcn;   //!
   TBranch        *b_hsf__value;   //!
   TBranch        *b_hsh__rcent;   //!
   TBranch        *b_hsh__fcent;   //!
   TBranch        *b_hsh__radius;   //!
   TBranch        *b_hsh__lambda;   //!
   TBranch        *b_hsh__fz0;   //!
   TBranch        *b_hst0;   //!
   TBranch        *b_hst0err;   //!
   TBranch        *b_hsn;   //!
   TBranch        *b_hsna;   //!
   TBranch        *b_ksf__value;   //!
   TBranch        *b_ksh__pars;   //!
   TBranch        *b_kst0;   //!
   TBranch        *b_kst0err;   //!
   TBranch        *b_ksm;   //!
   TBranch        *b_ksmerr;   //!
   TBranch        *b_ksn;   //!
   TBranch        *b_ksna;   //!
   TBranch        *b_kff__value;   //!
   TBranch        *b_kfh__pars;   //!
   TBranch        *b_kft0;   //!
   TBranch        *b_kft0err;   //!
   TBranch        *b_kfm;   //!
   TBranch        *b_kfmerr;   //!
   TBranch        *b_kfn;   //!
   TBranch        *b_kfna;   //!
   TBranch        *b_beamwt;   //!
   TBranch        *b_ndigitot;   //!
   TBranch        *b_mcgenmom;   //!
   TBranch        *b_mcentmom;   //!
   TBranch        *b_mcentpz;   //!
   TBranch        *b_mcentt0;   //!
   TBranch        *b_mcmidmom;   //!
   TBranch        *b_mcmidpz;   //!
   TBranch        *b_mcmidt0;   //!
   TBranch        *b_mcxitmom;   //!
   TBranch        *b_mcxitpz;   //!
   TBranch        *b_mcxitt0;   //!
   TBranch        *b_mch__rcent;   //!
   TBranch        *b_mch__fcent;   //!
   TBranch        *b_mch__radius;   //!
   TBranch        *b_mch__lambda;   //!
   TBranch        *b_mch__fz0;   //!
   TBranch        *b_pdg;   //!
   TBranch        *b_gen;   //!
   TBranch        *b_proc;   //!
   TBranch        *b_ndigi;   //!
   TBranch        *b_nkfprimary;   //!
   TBranch        *b_nksprimary;   //!
   TBranch        *b_nhsprimary;   //!
   TBranch        *b_ntcprimary;   //!
   TBranch        *b_kfnp;   //!
   TBranch        *b_kfnap;   //!
   TBranch        *b_ksnp;   //!
   TBranch        *b_ksnap;   //!
   TBranch        *b_hsnp;   //!
   TBranch        *b_hsnap;   //!
   TBranch        *b_tcnp;   //!

   TrkRecoDiag(TTree *tree,double norm);
   virtual ~TrkRecoDiag();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
// my variables
// this is copied from TrkFitFlag.hh, It should be #included FIXME!
  enum bit_type {hitsOK=0,circleOK,phizOK,helixOK,seedOK,kalmanOK,circleInit,phizInit,
    circleConverged,phizConverged,helixConverged,seedConverged,kalmanConverged,ntffbits};
  std::vector<unsigned> _tffval;
  double _norm; // normalization
  std::string _prefix; // prefix for root objects to avoid clash
  TH1F *_eff, *_acc; // signal efficiency , acceptance
// trigger study distributions
  TH1F *_hn, *_hna, *_hd0, *_hrmax, *_hrad, *_hlam, *_hmom;
  TH1F *_shn, *_shna, *_shd0, *_shrmax, *_shrad, *_shlam, *_shmom, *_ssmom, *_ssna;
  TH1F *_fhn, *_fhna, *_fhd0, *_fhrmax, *_fhrad, *_fhlam, *_fhmom, *_fsmom, *_fsna;

  TCanvas* _effcan;

  void createHistos();
  const char* title(const char*);
};

#endif

#ifdef TrkRecoDiag_cxx
TrkRecoDiag::~TrkRecoDiag()
{
   if (!_chain) return;
   delete _chain->GetCurrentFile();
}

Int_t TrkRecoDiag::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!_chain) return 0;
   return _chain->GetEntry(entry);
}
Long64_t TrkRecoDiag::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!_chain) return -5;
   Long64_t centry = _chain->LoadTree(entry);
   if (centry < 0) return centry;
   if (_chain->GetTreeNumber() != fCurrent) {
      fCurrent = _chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void TrkRecoDiag::Init(TTree *tree)
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
   _chain = tree;
   fCurrent = -1;
   _chain->SetMakeClass(1);

   _chain->SetBranchAddress("iev", &iev, &b_iev);
   _chain->SetBranchAddress("tct0", &tct0, &b_tct0);
   _chain->SetBranchAddress("tct0err", &tct0err, &b_tct0err);
   _chain->SetBranchAddress("tcn", &tcn, &b_tcn);
   _chain->SetBranchAddress("hsf._value", &hsf__value, &b_hsf__value);
   _chain->SetBranchAddress("hsh._rcent", &hsh__rcent, &b_hsh__rcent);
   _chain->SetBranchAddress("hsh._fcent", &hsh__fcent, &b_hsh__fcent);
   _chain->SetBranchAddress("hsh._radius", &hsh__radius, &b_hsh__radius);
   _chain->SetBranchAddress("hsh._lambda", &hsh__lambda, &b_hsh__lambda);
   _chain->SetBranchAddress("hsh._fz0", &hsh__fz0, &b_hsh__fz0);
   _chain->SetBranchAddress("hst0", &hst0, &b_hst0);
   _chain->SetBranchAddress("hst0err", &hst0err, &b_hst0err);
   _chain->SetBranchAddress("hsn", &hsn, &b_hsn);
   _chain->SetBranchAddress("hsna", &hsna, &b_hsna);
   _chain->SetBranchAddress("ksf._value", &ksf__value, &b_ksf__value);
   _chain->SetBranchAddress("ksh._pars[5]", ksh__pars, &b_ksh__pars);
   _chain->SetBranchAddress("kst0", &kst0, &b_kst0);
   _chain->SetBranchAddress("kst0err", &kst0err, &b_kst0err);
   _chain->SetBranchAddress("ksm", &ksm, &b_ksm);
   _chain->SetBranchAddress("ksmerr", &ksmerr, &b_ksmerr);
   _chain->SetBranchAddress("ksn", &ksn, &b_ksn);
   _chain->SetBranchAddress("ksna", &ksna, &b_ksna);
   _chain->SetBranchAddress("kff._value", &kff__value, &b_kff__value);
   _chain->SetBranchAddress("kfh._pars[5]", kfh__pars, &b_kfh__pars);
   _chain->SetBranchAddress("kft0", &kft0, &b_kft0);
   _chain->SetBranchAddress("kft0err", &kft0err, &b_kft0err);
   _chain->SetBranchAddress("kfm", &kfm, &b_kfm);
   _chain->SetBranchAddress("kfmerr", &kfmerr, &b_kfmerr);
   _chain->SetBranchAddress("kfn", &kfn, &b_kfn);
   _chain->SetBranchAddress("kfna", &kfna, &b_kfna);
   _chain->SetBranchAddress("beamwt", &beamwt, &b_beamwt);
   _chain->SetBranchAddress("ndigitot", &ndigitot, &b_ndigitot);
   _chain->SetBranchAddress("mcgenmom", &mcgenmom, &b_mcgenmom);
   _chain->SetBranchAddress("mcentmom", &mcentmom, &b_mcentmom);
   _chain->SetBranchAddress("mcentpz", &mcentpz, &b_mcentpz);
   _chain->SetBranchAddress("mcentt0", &mcentt0, &b_mcentt0);
   _chain->SetBranchAddress("mcmidmom", &mcmidmom, &b_mcmidmom);
   _chain->SetBranchAddress("mcmidpz", &mcmidpz, &b_mcmidpz);
   _chain->SetBranchAddress("mcmidt0", &mcmidt0, &b_mcmidt0);
   _chain->SetBranchAddress("mcxitmom", &mcxitmom, &b_mcxitmom);
   _chain->SetBranchAddress("mcxitpz", &mcxitpz, &b_mcxitpz);
   _chain->SetBranchAddress("mcxitt0", &mcxitt0, &b_mcxitt0);
   _chain->SetBranchAddress("mch._rcent", &mch__rcent, &b_mch__rcent);
   _chain->SetBranchAddress("mch._fcent", &mch__fcent, &b_mch__fcent);
   _chain->SetBranchAddress("mch._radius", &mch__radius, &b_mch__radius);
   _chain->SetBranchAddress("mch._lambda", &mch__lambda, &b_mch__lambda);
   _chain->SetBranchAddress("mch._fz0", &mch__fz0, &b_mch__fz0);
   _chain->SetBranchAddress("pdg", &pdg, &b_pdg);
   _chain->SetBranchAddress("gen", &gen, &b_gen);
   _chain->SetBranchAddress("proc", &proc, &b_proc);
   _chain->SetBranchAddress("ndigi", &ndigi, &b_ndigi);
   _chain->SetBranchAddress("nkfprimary", &nkfprimary, &b_nkfprimary);
   _chain->SetBranchAddress("nksprimary", &nksprimary, &b_nksprimary);
   _chain->SetBranchAddress("nhsprimary", &nhsprimary, &b_nhsprimary);
   _chain->SetBranchAddress("ntcprimary", &ntcprimary, &b_ntcprimary);
   _chain->SetBranchAddress("kfnp", &kfnp, &b_kfnp);
   _chain->SetBranchAddress("kfnap", &kfnap, &b_kfnap);
   _chain->SetBranchAddress("ksnp", &ksnp, &b_ksnp);
   _chain->SetBranchAddress("ksnap", &ksnap, &b_ksnap);
   _chain->SetBranchAddress("hsnp", &hsnp, &b_hsnp);
   _chain->SetBranchAddress("hsnap", &hsnap, &b_hsnap);
   _chain->SetBranchAddress("tcnp", &tcnp, &b_tcnp);
   Notify();
}

Bool_t TrkRecoDiag::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void TrkRecoDiag::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!_chain) return;
   _chain->Show(entry);
}
#endif // #ifdef TrkRecoDiag_cxx
