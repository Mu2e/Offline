//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon May 14 11:17:07 2018 by ROOT version 6.10/06
// from TChain TrkAna/trkana/
//////////////////////////////////////////////////////////

#ifndef HitAmbigQual_h
#define HitAmbigQual_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>

// Header file for the classes stored in the TTree if any.
#include "Math/GenVector/Cartesian3D.h"
#include "vector"
//#include "TrkDiag/inc/TrkStrawHitInfo.hh"
//#include "TrkDiag/inc/TrkStrawMatInfo.hh"
//#include "TrkDiag/inc/TrkStrawHitInfoMC.hh"

class HitAmbigQual {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
   static constexpr Int_t kMaxdemtsh = 94;
   static constexpr Int_t kMaxdemtsm = 154;
   static constexpr Int_t kMaxdemtshmc = 94;

   // Declaration of leaf types
   Int_t           evtinfo__eventid;
   Int_t           evtinfo__runid;
   Int_t           evtinfo__subrunid;
   Float_t         evtinfo__evtwt;
   Float_t         evtinfo__beamwt;
   Float_t         evtinfo__genwt;
   Int_t           evtinfo__nprotons;
   Int_t           hcnt__nsh;
   Int_t           hcnt__nesel;
   Int_t           hcnt__nrsel;
   Int_t           hcnt__ntsel;
   Int_t           hcnt__nbkg;
   Int_t           hcnt__nster;
   Int_t           hcnt__ntdiv;
   Int_t           hcnt__ntpk;
   Int_t           hcnt__nxt;
   Int_t           tcnt__ndem;
   Int_t           tcnt__nuem;
   Int_t           tcnt__ndmm;
   Int_t           tcnt__ndemc;
   Int_t           tcnt__ndemo;
   Int_t           tcnt__ndmmo;
   Int_t           dem__status;
   Int_t           dem__pdg;
   Int_t           dem__nhits;
   Int_t           dem__ndof;
   Int_t           dem__nactive;
   Int_t           dem__ndouble;
   Int_t           dem__ndactive;
   Int_t           dem__nnullambig;
   Int_t           dem__nmat;
   Int_t           dem__nmatactive;
   Int_t           dem__nbend;
   Float_t         dem__t0;
   Float_t         dem__t0err;
   Float_t         dem__chisq;
   Float_t         dem__con;
   Float_t         dem__radlen;
   Float_t         dem__firstflt;
   Float_t         dem__lastflt;
   Float_t         dem__startvalid;
   Float_t         dem__endvalid;
   Float_t         dem__trkqual;
   Float_t         dem__mom;
   Float_t         dem__momerr;
   Float_t         dem__fltlen;
   Float_t         dem__d0;
   Float_t         dem__p0;
   Float_t         dem__om;
   Float_t         dem__z0;
   Float_t         dem__td;
   Float_t         dem__d0err;
   Float_t         dem__p0err;
   Float_t         dem__omerr;
   Float_t         dem__z0err;
   Float_t         dem__tderr;
   Int_t           demtsh_;
   Bool_t          demtsh__active[kMaxdemtsh];   //[demtsh_]
   Bool_t          demtsh__dhit[kMaxdemtsh];   //[demtsh_]
   Bool_t          demtsh__dactive[kMaxdemtsh];   //[demtsh_]
   Int_t           demtsh__plane[kMaxdemtsh];   //[demtsh_]
   Int_t           demtsh__panel[kMaxdemtsh];   //[demtsh_]
   Int_t           demtsh__layer[kMaxdemtsh];   //[demtsh_]
   Int_t           demtsh__straw[kMaxdemtsh];   //[demtsh_]
   Int_t           demtsh__ambig[kMaxdemtsh];   //[demtsh_]
   Int_t           demtsh__driftend[kMaxdemtsh];   //[demtsh_]
   Float_t         demtsh__tdrift[kMaxdemtsh];   //[demtsh_]
   Float_t         demtsh__poca_fCoordinates_fX[kMaxdemtsh];   //[demtsh_]
   Float_t         demtsh__poca_fCoordinates_fY[kMaxdemtsh];   //[demtsh_]
   Float_t         demtsh__poca_fCoordinates_fZ[kMaxdemtsh];   //[demtsh_]
   Float_t         demtsh__resid[kMaxdemtsh];   //[demtsh_]
   Float_t         demtsh__residerr[kMaxdemtsh];   //[demtsh_]
   Float_t         demtsh__rdrift[kMaxdemtsh];   //[demtsh_]
   Float_t         demtsh__rdrifterr[kMaxdemtsh];   //[demtsh_]
   Float_t         demtsh__wdist[kMaxdemtsh];   //[demtsh_]
   Float_t         demtsh__werr[kMaxdemtsh];   //[demtsh_]
   Float_t         demtsh__trklen[kMaxdemtsh];   //[demtsh_]
   Float_t         demtsh__doca[kMaxdemtsh];   //[demtsh_]
   Float_t         demtsh__exerr[kMaxdemtsh];   //[demtsh_]
   Float_t         demtsh__penerr[kMaxdemtsh];   //[demtsh_]
   Float_t         demtsh__t0[kMaxdemtsh];   //[demtsh_]
   Float_t         demtsh__t0err[kMaxdemtsh];   //[demtsh_]
   Float_t         demtsh__ht[kMaxdemtsh];   //[demtsh_]
   Float_t         demtsh__hlen[kMaxdemtsh];   //[demtsh_]
   Float_t         demtsh__wdot[kMaxdemtsh];   //[demtsh_]
   Float_t         demtsh__bdot[kMaxdemtsh];   //[demtsh_]
   Float_t         demtsh__edep[kMaxdemtsh];   //[demtsh_]
   Float_t         demtsh__dx[kMaxdemtsh];   //[demtsh_]
   Int_t           demtsm_;
   Bool_t          demtsm__active[kMaxdemtsm];   //[demtsm_]
   Bool_t          demtsm__thit[kMaxdemtsm];   //[demtsm_]
   Bool_t          demtsm__thita[kMaxdemtsm];   //[demtsm_]
   Int_t           demtsm__plane[kMaxdemtsm];   //[demtsm_]
   Int_t           demtsm__panel[kMaxdemtsm];   //[demtsm_]
   Int_t           demtsm__layer[kMaxdemtsm];   //[demtsm_]
   Int_t           demtsm__straw[kMaxdemtsm];   //[demtsm_]
   Float_t         demtsm__doca[kMaxdemtsm];   //[demtsm_]
   Float_t         demtsm__tlen[kMaxdemtsm];   //[demtsm_]
   Float_t         demtsm__dp[kMaxdemtsm];   //[demtsm_]
   Float_t         demtsm__radlen[kMaxdemtsm];   //[demtsm_]
   Float_t         demtsm__sigMS[kMaxdemtsm];   //[demtsm_]
   Int_t           uem__status;
   Int_t           uem__pdg;
   Int_t           uem__nhits;
   Int_t           uem__ndof;
   Int_t           uem__nactive;
   Int_t           uem__ndouble;
   Int_t           uem__ndactive;
   Int_t           uem__nnullambig;
   Int_t           uem__nmat;
   Int_t           uem__nmatactive;
   Int_t           uem__nbend;
   Float_t         uem__t0;
   Float_t         uem__t0err;
   Float_t         uem__chisq;
   Float_t         uem__con;
   Float_t         uem__radlen;
   Float_t         uem__firstflt;
   Float_t         uem__lastflt;
   Float_t         uem__startvalid;
   Float_t         uem__endvalid;
   Float_t         uem__trkqual;
   Float_t         uem__mom;
   Float_t         uem__momerr;
   Float_t         uem__fltlen;
   Float_t         uem__d0;
   Float_t         uem__p0;
   Float_t         uem__om;
   Float_t         uem__z0;
   Float_t         uem__td;
   Float_t         uem__d0err;
   Float_t         uem__p0err;
   Float_t         uem__omerr;
   Float_t         uem__z0err;
   Float_t         uem__tderr;
   Int_t           dmm__status;
   Int_t           dmm__pdg;
   Int_t           dmm__nhits;
   Int_t           dmm__ndof;
   Int_t           dmm__nactive;
   Int_t           dmm__ndouble;
   Int_t           dmm__ndactive;
   Int_t           dmm__nnullambig;
   Int_t           dmm__nmat;
   Int_t           dmm__nmatactive;
   Int_t           dmm__nbend;
   Float_t         dmm__t0;
   Float_t         dmm__t0err;
   Float_t         dmm__chisq;
   Float_t         dmm__con;
   Float_t         dmm__radlen;
   Float_t         dmm__firstflt;
   Float_t         dmm__lastflt;
   Float_t         dmm__startvalid;
   Float_t         dmm__endvalid;
   Float_t         dmm__trkqual;
   Float_t         dmm__mom;
   Float_t         dmm__momerr;
   Float_t         dmm__fltlen;
   Float_t         dmm__d0;
   Float_t         dmm__p0;
   Float_t         dmm__om;
   Float_t         dmm__z0;
   Float_t         dmm__td;
   Float_t         dmm__d0err;
   Float_t         dmm__p0err;
   Float_t         dmm__omerr;
   Float_t         dmm__z0err;
   Float_t         dmm__tderr;
   Float_t         demc__dt;
   Float_t         demc__du;
   Float_t         demc__dv;
   Float_t         demc__ds;
   Float_t         demc__ep;
   Float_t         demc__uvchisq;
   Float_t         demc__tchisq;
   Float_t         demc__dtllr;
   Float_t         demc__epllr;
   Float_t         demc__eclust;
   Float_t         demc__tclust;
   Int_t           demc__section;
   Float_t         demc__cposx;
   Float_t         demc__cposy;
   Float_t         demc__cposz;
   Float_t         demc__tposx;
   Float_t         demc__tposy;
   Float_t         demc__tposz;
   Float_t         demc__tdirx;
   Float_t         demc__tdiry;
   Float_t         demc__tdirz;
   Float_t         demc__ttrk;
   Int_t           demmc_ndigi;
   Int_t           demmc_ndigigood;
   Int_t           demmc_nhits;
   Int_t           demmc_nactive;
   Int_t           demmc_ngood;
   Int_t           demmc_nambig;
   Int_t           demmc_pdg;
   Int_t           demmc_gen;
   Int_t           demmc_proc;
   Int_t           demmc_ppdg;
   Int_t           demmc_pgen;
   Int_t           demmc_pproc;
   Float_t         demmc_pmom;
   Float_t         demmcgen_t0;
   Float_t         demmcgen_mom;
   Float_t         demmcgen_x;
   Float_t         demmcgen_y;
   Float_t         demmcgen_z;
   Float_t         demmcgen_d0;
   Float_t         demmcgen_p0;
   Float_t         demmcgen_om;
   Float_t         demmcgen_z0;
   Float_t         demmcgen_td;
   Float_t         demmcent_t0;
   Float_t         demmcent_mom;
   Float_t         demmcent_x;
   Float_t         demmcent_y;
   Float_t         demmcent_z;
   Float_t         demmcent_d0;
   Float_t         demmcent_p0;
   Float_t         demmcent_om;
   Float_t         demmcent_z0;
   Float_t         demmcent_td;
   Float_t         demmcmid_t0;
   Float_t         demmcmid_mom;
   Float_t         demmcmid_x;
   Float_t         demmcmid_y;
   Float_t         demmcmid_z;
   Float_t         demmcmid_d0;
   Float_t         demmcmid_p0;
   Float_t         demmcmid_om;
   Float_t         demmcmid_z0;
   Float_t         demmcmid_td;
   Float_t         demmcxit_t0;
   Float_t         demmcxit_mom;
   Float_t         demmcxit_x;
   Float_t         demmcxit_y;
   Float_t         demmcxit_z;
   Float_t         demmcxit_d0;
   Float_t         demmcxit_p0;
   Float_t         demmcxit_om;
   Float_t         demmcxit_z0;
   Float_t         demmcxit_td;
   Int_t           demtshmc_;
   Int_t           demtshmc__pdg[kMaxdemtshmc];   //[demtshmc_]
   Int_t           demtshmc__gen[kMaxdemtshmc];   //[demtshmc_]
   Int_t           demtshmc__proc[kMaxdemtshmc];   //[demtshmc_]
   Int_t           demtshmc__rel[kMaxdemtshmc];   //[demtshmc_]
   Float_t         demtshmc__t0[kMaxdemtshmc];   //[demtshmc_]
   Float_t         demtshmc__ht[kMaxdemtshmc];   //[demtshmc_]
   Float_t         demtshmc__dist[kMaxdemtshmc];   //[demtshmc_]
   Float_t         demtshmc__len[kMaxdemtshmc];   //[demtshmc_]
   Float_t         demtshmc__edep[kMaxdemtshmc];   //[demtshmc_]
   Float_t         demtshmc__mom[kMaxdemtshmc];   //[demtshmc_]
   Float_t         demtshmc__r[kMaxdemtshmc];   //[demtshmc_]
   Float_t         demtshmc__phi[kMaxdemtshmc];   //[demtshmc_]
   Int_t           demtshmc__ambig[kMaxdemtshmc];   //[demtshmc_]
   Bool_t          demtshmc__xtalk[kMaxdemtshmc];   //[demtshmc_]

   // List of branches
   TBranch        *b_evtinfo_;   //!
   TBranch        *b_hcnt_;   //!
   TBranch        *b_tcnt_;   //!
   TBranch        *b_dem_;   //!
   TBranch        *b_demtsh_;   //!
   TBranch        *b_demtsh__active;   //!
   TBranch        *b_demtsh__dhit;   //!
   TBranch        *b_demtsh__dactive;   //!
   TBranch        *b_demtsh__plane;   //!
   TBranch        *b_demtsh__panel;   //!
   TBranch        *b_demtsh__layer;   //!
   TBranch        *b_demtsh__straw;   //!
   TBranch        *b_demtsh__ambig;   //!
   TBranch        *b_demtsh__driftend;   //!
   TBranch        *b_demtsh__tdrift;   //!
   TBranch        *b_demtsh__poca_fCoordinates_fX;   //!
   TBranch        *b_demtsh__poca_fCoordinates_fY;   //!
   TBranch        *b_demtsh__poca_fCoordinates_fZ;   //!
   TBranch        *b_demtsh__resid;   //!
   TBranch        *b_demtsh__residerr;   //!
   TBranch        *b_demtsh__rdrift;   //!
   TBranch        *b_demtsh__rdrifterr;   //!
   TBranch        *b_demtsh__wdist;   //!
   TBranch        *b_demtsh__werr;   //!
   TBranch        *b_demtsh__trklen;   //!
   TBranch        *b_demtsh__doca;   //!
   TBranch        *b_demtsh__exerr;   //!
   TBranch        *b_demtsh__penerr;   //!
   TBranch        *b_demtsh__t0;   //!
   TBranch        *b_demtsh__t0err;   //!
   TBranch        *b_demtsh__ht;   //!
   TBranch        *b_demtsh__hlen;   //!
   TBranch        *b_demtsh__wdot;   //!
   TBranch        *b_demtsh__bdot;   //!
   TBranch        *b_demtsh__edep;   //!
   TBranch        *b_demtsh__dx;   //!
   TBranch        *b_demtsm_;   //!
   TBranch        *b_demtsm__active;   //!
   TBranch        *b_demtsm__thit;   //!
   TBranch        *b_demtsm__thita;   //!
   TBranch        *b_demtsm__plane;   //!
   TBranch        *b_demtsm__panel;   //!
   TBranch        *b_demtsm__layer;   //!
   TBranch        *b_demtsm__straw;   //!
   TBranch        *b_demtsm__doca;   //!
   TBranch        *b_demtsm__tlen;   //!
   TBranch        *b_demtsm__dp;   //!
   TBranch        *b_demtsm__radlen;   //!
   TBranch        *b_demtsm__sigMS;   //!
   TBranch        *b_uem_;   //!
   TBranch        *b_dmm_;   //!
   TBranch        *b_demc_;   //!
   TBranch        *b_demmc;   //!
   TBranch        *b_demmcgen;   //!
   TBranch        *b_demmcent;   //!
   TBranch        *b_demmcmid;   //!
   TBranch        *b_demmcxit;   //!
   TBranch        *b_demtshmc_;   //!
   TBranch        *b_demtshmc__pdg;   //!
   TBranch        *b_demtshmc__gen;   //!
   TBranch        *b_demtshmc__proc;   //!
   TBranch        *b_demtshmc__rel;   //!
   TBranch        *b_demtshmc__t0;   //!
   TBranch        *b_demtshmc__ht;   //!
   TBranch        *b_demtshmc__dist;   //!
   TBranch        *b_demtshmc__len;   //!
   TBranch        *b_demtshmc__edep;   //!
   TBranch        *b_demtshmc__mom;   //!
   TBranch        *b_demtshmc__r;   //!
   TBranch        *b_demtshmc__phi;   //!
   TBranch        *b_demtshmc__ambig;   //!
   TBranch        *b_demtshmc__xtalk;   //!

   TH1F *_paplane, *_papanel, *_pastation;
   TH1F *_spaplane, *_spapanel, *_spastation;
   TH1F *_spadplane, *_spadpanel, *_spadstation;
   TH1F *_spasplane, *_spaspanel, *_spasstation;
   TH1F *_spasdplane, *_spasdpanel, *_spasdstation;
   TH1F *_spasfplane, *_spasfpanel, *_spasfstation;
   TH1F *_spasdfplane, *_spasdfpanel, *_spasdfstation;
   TH1F *_pdplane, *_pdpanel, *_pdstation;

   TH2F *_dmplane, *_dmpanel, *_dmstation;
   TH2F *_dmtplane, *_dmtpanel, *_dmtstation;


   HitAmbigQual(TTree *tree=0);
   virtual ~HitAmbigQual();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef HitAmbigQual_cxx
HitAmbigQual::HitAmbigQual(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/data/7852374/00/00000/nts.brownd.TAM.TAM.004001_00000000.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/data/7852374/00/00000/nts.brownd.TAM.TAM.004001_00000000.root");
      }
      f->GetObject("TrkAna/trkana",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("TrkAna/trkana","");
      chain->Add("/data/7852374/00/00000/nts.brownd.TAM.TAM.004001_00000000.root/TrkAna/trkana");
      chain->Add("/data/7852374/00/00001/nts.brownd.TAM.TAM.004001_00000001.root/TrkAna/trkana");
      chain->Add("/data/7852374/00/00002/nts.brownd.TAM.TAM.004001_00000002.root/TrkAna/trkana");
      chain->Add("/data/7852374/00/00003/nts.brownd.TAM.TAM.004001_00000003.root/TrkAna/trkana");
      chain->Add("/data/7852374/00/00004/nts.brownd.TAM.TAM.004001_00000004.root/TrkAna/trkana");
      chain->Add("/data/7852374/00/00005/nts.brownd.TAM.TAM.004001_00000005.root/TrkAna/trkana");
      chain->Add("/data/7852374/00/00006/nts.brownd.TAM.TAM.004001_00000006.root/TrkAna/trkana");
      chain->Add("/data/7852374/00/00007/nts.brownd.TAM.TAM.004001_00000007.root/TrkAna/trkana");
      chain->Add("/data/7852374/00/00008/nts.brownd.TAM.TAM.004001_00000008.root/TrkAna/trkana");
      chain->Add("/data/7852374/00/00009/nts.brownd.TAM.TAM.004001_00000009.root/TrkAna/trkana");
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

HitAmbigQual::~HitAmbigQual()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t HitAmbigQual::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t HitAmbigQual::LoadTree(Long64_t entry)
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

void HitAmbigQual::Init(TTree *tree)
{

  _papanel = new TH1F("papanel","Pair Ambiguity (panel)",100,-7.0,7.0);
  _paplane = new TH1F("paplane","Pair Ambiguity (plane)",100,-7.0,7.0);
  _pastation = new TH1F("pastation","Pair Ambiguity (station)",100,-7.0,7.0);

  _spaspanel = new TH1F("spaspanel","Signed Sum Pair Ambiguity (panel)",100,0.0,100.0);
  _spasplane = new TH1F("spasplane","Signed Sum Pair Ambiguity (plane)",100,0.0,200.0);
  _spasstation = new TH1F("spasstation","Signed Sum Pair Ambiguity (station)",100,0.0,400.0);

  _spasdpanel = new TH1F("spasdpanel","DNorm Signed Sum Pair Ambiguity (panel)",100,0.0,10.0);
  _spasdplane = new TH1F("spasdplane","DNorm Signed Sum Pair Ambiguity (plane)",100,0.0,20.0);
  _spasdstation = new TH1F("spasdstation","DNorm Signed Sum Pair Ambiguity (station)",100,0.0,40.0);

  _spasfpanel = new TH1F("spasfpanel","Signed Sum Pair Ambiguity Fraction (panel)",100,0.0,1.0);
  _spasfplane = new TH1F("spasfplane","Signed Sum Pair Ambiguity Fraction (plane)",100,0.0,1.0);
  _spasfstation = new TH1F("spasfstation","Signed Sum Pair Ambiguity Fraction (station)",100,0.0,1.0);

  _spasdfpanel = new TH1F("spasdfpanel","DNorm Signed Sum Pair Ambiguity Fraction (panel)",100,0.0,1.0);
  _spasdfplane = new TH1F("spasdfplane","DNorm Signed Sum Pair Ambiguity Fraction (plane)",100,0.0,1.0);
  _spasdfstation = new TH1F("spasdfstation","DNorm Signed Sum Pair Ambiguity Fraction (station)",100,0.0,1.0);

  _spapanel = new TH1F("spapanel","Sum Pair Ambiguity (panel)",100,0.0,100.0);
  _spaplane = new TH1F("spaplane","Sum Pair Ambiguity (plane)",100,0.0,200.0);
  _spastation = new TH1F("spastation","Sum Pair Ambiguity (station)",100,0.0,400.0);

  _spadpanel = new TH1F("spadpanel","DNorm Sum Pair Ambiguity (panel)",100,0.0,10.0);
  _spadplane = new TH1F("spadplane","DNorm Sum Pair Ambiguity (plane)",100,0.0,20.0);
  _spadstation = new TH1F("spadstation","DNorm Sum Pair Ambiguity (station)",100,0.0,40.0);

  _pdpanel = new TH1F("pdpanel","Pair Distance (panel)",100,0.0,20.0);
  _pdplane = new TH1F("pdplane","Pair Distance (plane)",100,0.0,70.0);
  _pdstation = new TH1F("pdstation","Pair Distance (station)",100,0.0,200.0);

  _dmpanel = new TH2F("dmpanel","Mom Diff vs DNorm Ambig Frac (panel)",50,0.0,1.2,50,0.0,5.0);
  _dmplane = new TH2F("dmplane","Mom Diff vs DNorm Ambig Frac (plane)",50,0.0,1.2,50,0.0,5.0);
  _dmstation = new TH2F("dmstation","Mom Diff vs DNorm Ambig Frac (station)",50,0.0,1.2,50,0.0,5.0);

  _dmtpanel = new TH2F("dmtpanel","Mom Diff vs DNorm Ambig Frac Trkqual (panel)",50,0.0,1.2,50,0.0,5.0);
  _dmtplane = new TH2F("dmtplane","Mom Diff vs DNorm Ambig Frac Trkqual (plane)",50,0.0,1.2,50,0.0,5.0);
  _dmtstation = new TH2F("dmtstation","Mom Diff vs DNorm Ambig Frac Trkqual (station)",50,0.0,1.2,50,0.0,5.0);

 
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

   fChain->SetBranchAddress("evtinfo.", &evtinfo__eventid, &b_evtinfo_);
   fChain->SetBranchAddress("hcnt.", &hcnt__nsh, &b_hcnt_);
   fChain->SetBranchAddress("tcnt.", &tcnt__ndem, &b_tcnt_);
   fChain->SetBranchAddress("dem.", &dem__status, &b_dem_);
   fChain->SetBranchAddress("demtsh", &demtsh_, &b_demtsh_);
   fChain->SetBranchAddress("demtsh._active", demtsh__active, &b_demtsh__active);
   fChain->SetBranchAddress("demtsh._dhit", demtsh__dhit, &b_demtsh__dhit);
   fChain->SetBranchAddress("demtsh._dactive", demtsh__dactive, &b_demtsh__dactive);
   fChain->SetBranchAddress("demtsh._plane", demtsh__plane, &b_demtsh__plane);
   fChain->SetBranchAddress("demtsh._panel", demtsh__panel, &b_demtsh__panel);
   fChain->SetBranchAddress("demtsh._layer", demtsh__layer, &b_demtsh__layer);
   fChain->SetBranchAddress("demtsh._straw", demtsh__straw, &b_demtsh__straw);
   fChain->SetBranchAddress("demtsh._ambig", demtsh__ambig, &b_demtsh__ambig);
   fChain->SetBranchAddress("demtsh._driftend", demtsh__driftend, &b_demtsh__driftend);
   fChain->SetBranchAddress("demtsh._tdrift", demtsh__tdrift, &b_demtsh__tdrift);
   fChain->SetBranchAddress("demtsh._poca.fCoordinates.fX", demtsh__poca_fCoordinates_fX, &b_demtsh__poca_fCoordinates_fX);
   fChain->SetBranchAddress("demtsh._poca.fCoordinates.fY", demtsh__poca_fCoordinates_fY, &b_demtsh__poca_fCoordinates_fY);
   fChain->SetBranchAddress("demtsh._poca.fCoordinates.fZ", demtsh__poca_fCoordinates_fZ, &b_demtsh__poca_fCoordinates_fZ);
   fChain->SetBranchAddress("demtsh._resid", demtsh__resid, &b_demtsh__resid);
   fChain->SetBranchAddress("demtsh._residerr", demtsh__residerr, &b_demtsh__residerr);
   fChain->SetBranchAddress("demtsh._rdrift", demtsh__rdrift, &b_demtsh__rdrift);
   fChain->SetBranchAddress("demtsh._rdrifterr", demtsh__rdrifterr, &b_demtsh__rdrifterr);
   fChain->SetBranchAddress("demtsh._wdist", demtsh__wdist, &b_demtsh__wdist);
   fChain->SetBranchAddress("demtsh._werr", demtsh__werr, &b_demtsh__werr);
   fChain->SetBranchAddress("demtsh._trklen", demtsh__trklen, &b_demtsh__trklen);
   fChain->SetBranchAddress("demtsh._doca", demtsh__doca, &b_demtsh__doca);
   fChain->SetBranchAddress("demtsh._exerr", demtsh__exerr, &b_demtsh__exerr);
   fChain->SetBranchAddress("demtsh._penerr", demtsh__penerr, &b_demtsh__penerr);
   fChain->SetBranchAddress("demtsh._t0", demtsh__t0, &b_demtsh__t0);
   fChain->SetBranchAddress("demtsh._t0err", demtsh__t0err, &b_demtsh__t0err);
   fChain->SetBranchAddress("demtsh._ht", demtsh__ht, &b_demtsh__ht);
   fChain->SetBranchAddress("demtsh._hlen", demtsh__hlen, &b_demtsh__hlen);
   fChain->SetBranchAddress("demtsh._wdot", demtsh__wdot, &b_demtsh__wdot);
   fChain->SetBranchAddress("demtsh._bdot", demtsh__bdot, &b_demtsh__bdot);
   fChain->SetBranchAddress("demtsh._edep", demtsh__edep, &b_demtsh__edep);
   fChain->SetBranchAddress("demtsh._dx", demtsh__dx, &b_demtsh__dx);
   fChain->SetBranchAddress("demtsm", &demtsm_, &b_demtsm_);
   fChain->SetBranchAddress("demtsm._active", demtsm__active, &b_demtsm__active);
   fChain->SetBranchAddress("demtsm._thit", demtsm__thit, &b_demtsm__thit);
   fChain->SetBranchAddress("demtsm._thita", demtsm__thita, &b_demtsm__thita);
   fChain->SetBranchAddress("demtsm._plane", demtsm__plane, &b_demtsm__plane);
   fChain->SetBranchAddress("demtsm._panel", demtsm__panel, &b_demtsm__panel);
   fChain->SetBranchAddress("demtsm._layer", demtsm__layer, &b_demtsm__layer);
   fChain->SetBranchAddress("demtsm._straw", demtsm__straw, &b_demtsm__straw);
   fChain->SetBranchAddress("demtsm._doca", demtsm__doca, &b_demtsm__doca);
   fChain->SetBranchAddress("demtsm._tlen", demtsm__tlen, &b_demtsm__tlen);
   fChain->SetBranchAddress("demtsm._dp", demtsm__dp, &b_demtsm__dp);
   fChain->SetBranchAddress("demtsm._radlen", demtsm__radlen, &b_demtsm__radlen);
   fChain->SetBranchAddress("demtsm._sigMS", demtsm__sigMS, &b_demtsm__sigMS);
   fChain->SetBranchAddress("uem.", &uem__status, &b_uem_);
   fChain->SetBranchAddress("dmm.", &dmm__status, &b_dmm_);
   fChain->SetBranchAddress("demc.", &demc__dt, &b_demc_);
   fChain->SetBranchAddress("demmc", &demmc_ndigi, &b_demmc);
   fChain->SetBranchAddress("demmcgen", &demmcgen_t0, &b_demmcgen);
   fChain->SetBranchAddress("demmcent", &demmcent_t0, &b_demmcent);
   fChain->SetBranchAddress("demmcmid", &demmcmid_t0, &b_demmcmid);
   fChain->SetBranchAddress("demmcxit", &demmcxit_t0, &b_demmcxit);
   fChain->SetBranchAddress("demtshmc", &demtshmc_, &b_demtshmc_);
   fChain->SetBranchAddress("demtshmc._pdg", demtshmc__pdg, &b_demtshmc__pdg);
   fChain->SetBranchAddress("demtshmc._gen", demtshmc__gen, &b_demtshmc__gen);
   fChain->SetBranchAddress("demtshmc._proc", demtshmc__proc, &b_demtshmc__proc);
   fChain->SetBranchAddress("demtshmc._rel", demtshmc__rel, &b_demtshmc__rel);
   fChain->SetBranchAddress("demtshmc._t0", demtshmc__t0, &b_demtshmc__t0);
   fChain->SetBranchAddress("demtshmc._ht", demtshmc__ht, &b_demtshmc__ht);
   fChain->SetBranchAddress("demtshmc._dist", demtshmc__dist, &b_demtshmc__dist);
   fChain->SetBranchAddress("demtshmc._len", demtshmc__len, &b_demtshmc__len);
   fChain->SetBranchAddress("demtshmc._edep", demtshmc__edep, &b_demtshmc__edep);
   fChain->SetBranchAddress("demtshmc._mom", demtshmc__mom, &b_demtshmc__mom);
   fChain->SetBranchAddress("demtshmc._r", demtshmc__r, &b_demtshmc__r);
   fChain->SetBranchAddress("demtshmc._phi", demtshmc__phi, &b_demtshmc__phi);
   fChain->SetBranchAddress("demtshmc._ambig", demtshmc__ambig, &b_demtshmc__ambig);
   fChain->SetBranchAddress("demtshmc._xtalk", demtshmc__xtalk, &b_demtshmc__xtalk);
   Notify();
}

Bool_t HitAmbigQual::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void HitAmbigQual::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t HitAmbigQual::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef HitAmbigQual_cxx
