#ifndef RobustMultiHelixFinderTypes_hh
#define RobustMultiHelixFinderTypes_hh

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "TTree.h"

namespace art {class Event;};

namespace mu2e {

  class ComboHitCollection;
  namespace RobustMultiHelixFinderTypes {

    struct Config{
        fhicl::Atom<std::string>   tool_type{             fhicl::Name("tool_type"),             fhicl::Comment("Needed fot backward compatibility"), "" };
        fhicl::Atom<bool>          mcDiag{                fhicl::Name("MCDiag"),                fhicl::Comment("Switch to perform MC diag"), true };
        fhicl::Atom<art::InputTag> strawDigiMCCollection{ fhicl::Name("StrawDigiMCCollection"), fhicl::Comment("StrawDigiMC collection name"), "makeSD" };
    };


    struct Data_t {

        enum  {kMaxHits = 8192};

        const art::Event*         event_;
        const ComboHitCollection* chcol_;

        Int_t   iev_;
        Int_t   Nch_,chNhit_[kMaxHits],chPdg_[kMaxHits],chCrCode_[kMaxHits],chSimId_[kMaxHits],chUId_[kMaxHits];
        Float_t chTime_[kMaxHits],chPhi_[kMaxHits],chRad_[kMaxHits],chX_[kMaxHits],chY_[kMaxHits],chZ_[kMaxHits];
        Float_t chTerr_[kMaxHits],chWerr_[kMaxHits],chWDX_[kMaxHits],chWDY_[kMaxHits];

        Int_t   Nhel_,helhel_[128],helnhi_[128];
        Float_t helrad_[128],helrcen_[128],helfcen_[128],hellam_[128],helfz0_[128],helchi2_[128],heldzdt_[128];
        std::vector<std::vector<int>> helhits_;
        std::vector<std::vector<uint16_t>> chStrawId_;

        void  reset() {Nch_=Nhel_=0;chStrawId_.clear();helhits_.clear();}
    };
  }
}

#endif
