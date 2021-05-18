//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Apr 12 10:55:52 2013 by ROOT version 5.34/05
// from TTree tstep/Stepper tree
// found on file: g4studyCalo_01.root
//////////////////////////////////////////////////////////

#ifndef g4studyCalo01Selector_h
#define g4studyCalo01Selector_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class g4studyCalo01Selector : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   Int_t           run;
   Int_t           evt;
   Int_t           trk;
   Int_t           vol;
   Int_t           pdg;
   Float_t         time;
   Float_t         gtime;
   Float_t         x;
   Float_t         y;
   Float_t         z;
   Float_t         px;
   Float_t         py;
   Float_t         pz;
   Float_t         p;
   Float_t         ke;
   Float_t         tedep;
   Float_t         niedep;
   Int_t           endcode;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_evt;   //!
   TBranch        *b_trk;   //!
   TBranch        *b_vol;   //!
   TBranch        *b_pdg;   //!
   TBranch        *b_time;   //!
   TBranch        *b_gtime;   //!
   TBranch        *b_x;   //!
   TBranch        *b_y;   //!
   TBranch        *b_z;   //!
   TBranch        *b_px;   //!
   TBranch        *b_py;   //!
   TBranch        *b_pz;   //!
   TBranch        *b_p;   //!
   TBranch        *b_ke;   //!
   TBranch        *b_tedep;   //!
   TBranch        *b_niedep;   //!
   TBranch        *b_endcode;   //!

   g4studyCalo01Selector(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~g4studyCalo01Selector() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(g4studyCalo01Selector,0);
};

#endif

#ifdef g4studyCalo01Selector_cxx
void g4studyCalo01Selector::Init(TTree *tree)
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
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("evt", &evt, &b_evt);
   fChain->SetBranchAddress("trk", &trk, &b_trk);
   fChain->SetBranchAddress("vol", &vol, &b_vol);
   fChain->SetBranchAddress("pdg", &pdg, &b_pdg);
   fChain->SetBranchAddress("time", &time, &b_time);
   fChain->SetBranchAddress("gtime", &gtime, &b_gtime);
   fChain->SetBranchAddress("x", &x, &b_x);
   fChain->SetBranchAddress("y", &y, &b_y);
   fChain->SetBranchAddress("z", &z, &b_z);
   fChain->SetBranchAddress("px", &px, &b_px);
   fChain->SetBranchAddress("py", &py, &b_py);
   fChain->SetBranchAddress("pz", &pz, &b_pz);
   fChain->SetBranchAddress("p", &p, &b_p);
   fChain->SetBranchAddress("ke", &ke, &b_ke);
   fChain->SetBranchAddress("tedep", &tedep, &b_tedep);
   fChain->SetBranchAddress("niedep", &niedep, &b_niedep);
   fChain->SetBranchAddress("endcode", &endcode, &b_endcode);
}

Bool_t g4studyCalo01Selector::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef g4studyCalo01Selector_cxx
