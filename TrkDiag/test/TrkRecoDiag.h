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
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
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
   Float_t         kschisq;
   Float_t         kscon;
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
   Float_t         kfchisq;
   Float_t         kfcon;
   Float_t         kftq;
   Int_t           kfn;
   Int_t           kfna;
   Int_t           kfns;
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
   TBranch        *b_kschisq;   //!
   TBranch        *b_kscon;   //!
   TBranch        *b_ksn;   //!
   TBranch        *b_ksna;   //!
   TBranch        *b_kff__value;   //!
   TBranch        *b_kfh__pars;   //!
   TBranch        *b_kft0;   //!
   TBranch        *b_kft0err;   //!
   TBranch        *b_kfm;   //!
   TBranch        *b_kfmerr;   //!
   TBranch        *b_kfchisq;   //!
   TBranch        *b_kfcon;   //!
   TBranch        *b_kftq;   //!
   TBranch        *b_kfn;   //!
   TBranch        *b_kfna;   //!
   TBranch        *b_kfns;   //!
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
  TH1F *_shn, *_shna, *_shd0, *_shrmax, *_shrad, *_shlam, *_shmom;
  TH1F *_fhn, *_fhna, *_fhd0, *_fhrmax, *_fhrad, *_fhlam, *_fhmom;
  TH1F *_phn, *_phna, *_phd0, *_phrmax, *_phrad, *_phlam, *_phmom;
  TH1F *_ssmom, *_ssna, *_sschisq, *_ssmerr;
  TH1F *_fsmom, *_fsna, *_fschisq, *_fsmerr; 
  TH1F *_psmom, *_psna, *_pschisq, *_psmerr; 
  
  TCanvas* _effcan;

  // cuts for defining 'trigger' selections
  int _tcseln;
  int _hseln;
  double _hselminm, _hselmaxm;
  double _sselminm, _sselmaxm, _sselmerr, _sselchi2;
  double _pseltq, _pselminm, _pselmaxm;
  double _mincost, _maxcost;
  double _mint0;
  double _tdlow, _tdhigh;


  void createHistos();
  const char* title(const char*);
};

#endif

#ifdef TrkRecoDiag_cxx
TrkRecoDiag::~TrkRecoDiag()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t TrkRecoDiag::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t TrkRecoDiag::LoadTree(Long64_t entry)
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
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("iev", &iev, &b_iev);
   fChain->SetBranchAddress("tct0", &tct0, &b_tct0);
   fChain->SetBranchAddress("tct0err", &tct0err, &b_tct0err);
   fChain->SetBranchAddress("tcn", &tcn, &b_tcn);
   fChain->SetBranchAddress("hsf._value", &hsf__value, &b_hsf__value);
   fChain->SetBranchAddress("hsh._rcent", &hsh__rcent, &b_hsh__rcent);
   fChain->SetBranchAddress("hsh._fcent", &hsh__fcent, &b_hsh__fcent);
   fChain->SetBranchAddress("hsh._radius", &hsh__radius, &b_hsh__radius);
   fChain->SetBranchAddress("hsh._lambda", &hsh__lambda, &b_hsh__lambda);
   fChain->SetBranchAddress("hsh._fz0", &hsh__fz0, &b_hsh__fz0);
   fChain->SetBranchAddress("hst0", &hst0, &b_hst0);
   fChain->SetBranchAddress("hst0err", &hst0err, &b_hst0err);
   fChain->SetBranchAddress("hsn", &hsn, &b_hsn);
   fChain->SetBranchAddress("hsna", &hsna, &b_hsna);
   fChain->SetBranchAddress("ksf._value", &ksf__value, &b_ksf__value);
   fChain->SetBranchAddress("ksh._pars[5]", ksh__pars, &b_ksh__pars);
   fChain->SetBranchAddress("kst0", &kst0, &b_kst0);
   fChain->SetBranchAddress("kst0err", &kst0err, &b_kst0err);
   fChain->SetBranchAddress("ksm", &ksm, &b_ksm);
   fChain->SetBranchAddress("ksmerr", &ksmerr, &b_ksmerr);
   fChain->SetBranchAddress("kschisq", &kschisq, &b_kschisq);
   fChain->SetBranchAddress("kscon", &kscon, &b_kscon);
   fChain->SetBranchAddress("ksn", &ksn, &b_ksn);
   fChain->SetBranchAddress("ksna", &ksna, &b_ksna);
   fChain->SetBranchAddress("kff._value", &kff__value, &b_kff__value);
   fChain->SetBranchAddress("kfh._pars[5]", kfh__pars, &b_kfh__pars);
   fChain->SetBranchAddress("kft0", &kft0, &b_kft0);
   fChain->SetBranchAddress("kft0err", &kft0err, &b_kft0err);
   fChain->SetBranchAddress("kfm", &kfm, &b_kfm);
   fChain->SetBranchAddress("kfmerr", &kfmerr, &b_kfmerr);
   fChain->SetBranchAddress("kfchisq", &kfchisq, &b_kfchisq);
   fChain->SetBranchAddress("kfcon", &kfcon, &b_kfcon);
   fChain->SetBranchAddress("kftq", &kftq, &b_kftq);
   fChain->SetBranchAddress("kfn", &kfn, &b_kfn);
   fChain->SetBranchAddress("kfna", &kfna, &b_kfna);
   fChain->SetBranchAddress("beamwt", &beamwt, &b_beamwt);
   fChain->SetBranchAddress("ndigitot", &ndigitot, &b_ndigitot);
   fChain->SetBranchAddress("mcgenmom", &mcgenmom, &b_mcgenmom);
   fChain->SetBranchAddress("mcentmom", &mcentmom, &b_mcentmom);
   fChain->SetBranchAddress("mcentpz", &mcentpz, &b_mcentpz);
   fChain->SetBranchAddress("mcentt0", &mcentt0, &b_mcentt0);
   fChain->SetBranchAddress("mcmidmom", &mcmidmom, &b_mcmidmom);
   fChain->SetBranchAddress("mcmidpz", &mcmidpz, &b_mcmidpz);
   fChain->SetBranchAddress("mcmidt0", &mcmidt0, &b_mcmidt0);
   fChain->SetBranchAddress("mcxitmom", &mcxitmom, &b_mcxitmom);
   fChain->SetBranchAddress("mcxitpz", &mcxitpz, &b_mcxitpz);
   fChain->SetBranchAddress("mcxitt0", &mcxitt0, &b_mcxitt0);
   fChain->SetBranchAddress("mch._rcent", &mch__rcent, &b_mch__rcent);
   fChain->SetBranchAddress("mch._fcent", &mch__fcent, &b_mch__fcent);
   fChain->SetBranchAddress("mch._radius", &mch__radius, &b_mch__radius);
   fChain->SetBranchAddress("mch._lambda", &mch__lambda, &b_mch__lambda);
   fChain->SetBranchAddress("mch._fz0", &mch__fz0, &b_mch__fz0);
   fChain->SetBranchAddress("pdg", &pdg, &b_pdg);
   fChain->SetBranchAddress("gen", &gen, &b_gen);
   fChain->SetBranchAddress("proc", &proc, &b_proc);
   fChain->SetBranchAddress("ndigi", &ndigi, &b_ndigi);
   fChain->SetBranchAddress("nkfprimary", &nkfprimary, &b_nkfprimary);
   fChain->SetBranchAddress("nksprimary", &nksprimary, &b_nksprimary);
   fChain->SetBranchAddress("nhsprimary", &nhsprimary, &b_nhsprimary);
   fChain->SetBranchAddress("ntcprimary", &ntcprimary, &b_ntcprimary);
   fChain->SetBranchAddress("kfnp", &kfnp, &b_kfnp);
   fChain->SetBranchAddress("kfnap", &kfnap, &b_kfnap);
   fChain->SetBranchAddress("ksnp", &ksnp, &b_ksnp);
   fChain->SetBranchAddress("ksnap", &ksnap, &b_ksnap);
   fChain->SetBranchAddress("hsnp", &hsnp, &b_hsnp);
   fChain->SetBranchAddress("hsnap", &hsnap, &b_hsnap);
   fChain->SetBranchAddress("tcnp", &tcnp, &b_tcnp);
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
   if (!fChain) return;
   fChain->Show(entry);
}
#endif // #ifdef TrkRecoDiag_cxx
