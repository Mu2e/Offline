#ifndef TimeAndPhiClusterFinderTypes_hh
#define TimeAndPhiClusterFinderTypes_hh

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "TTree.h"

namespace art {
  class Event;
};


namespace mu2e {

  class ComboHitCollection;
  
  namespace TimeAndPhiClusterFinderTypes {
  
    struct Config
    {
       fhicl::Atom<std::string>   tool_type{             fhicl::Name("tool_type"),             fhicl::Comment("Needed fot backward compatibility"), "" }; 
       fhicl::Atom<bool>          mcDiag{                fhicl::Name("MCDiag"),                fhicl::Comment("Switch to perform MC diag"), true }; 
       fhicl::Atom<art::InputTag> strawDigiMCCollection{ fhicl::Name("StrawDigiMCCollection"), fhicl::Comment("StrawDigiMC collection name"), "makeSD" };
    };       
    
    
    struct Data_t {
      
      enum  {kMaxHits = 8192};

      const art::Event*         event_;
      const ComboHitCollection* chcol_;

      Int_t   iev_;
      Int_t   Nch_,chSel_[kMaxHits];
      Int_t   Ncal_;
      Int_t   nhit1_,nclu1_[kMaxHits],hitIdx1_[kMaxHits];
      Int_t   nhit2_,nclu2_[kMaxHits],hitIdx2_[kMaxHits];
      Int_t   nhit3_,nclu3_[kMaxHits],hitIdx3_[kMaxHits];
      Float_t calTime_[kMaxHits], chTime_[kMaxHits], clu2Time_[kMaxHits],clu2Phi_[kMaxHits],clu2Rad_[kMaxHits];
      Int_t   nhitMVA_,icluMVA_, ncluMVA_[kMaxHits],hitIdxMVA_[kMaxHits],MVAsig_[kMaxHits],MVAsel_[kMaxHits];
      Float_t MVAvar1_[kMaxHits],MVAvar2_[kMaxHits],MVAvar3_[kMaxHits],MVAvar4_[kMaxHits],MVAvar5_[kMaxHits];      
    
      void  reset() {Nch_=nhit1_=nhit2_=nhit3_=nhitMVA_=icluMVA_=Ncal_=0;}
    };
  
  }
}

#endif
