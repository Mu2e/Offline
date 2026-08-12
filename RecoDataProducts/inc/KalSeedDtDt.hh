//
// Class representing the fit information of a dt_{straw hit time} / dt_{track hit time of closest approach}
// This is useful for PID, as the wrong hypothesis will typically have a slope
//  Original Author: Michael MacKenzie, 2026
//
#ifndef RecoDataProducts_KalSeedDtDt_HH
#define RecoDataProducts_KalSeedDtDt_HH
#include <Rtypes.h>
namespace mu2e {
  struct KalSeedDtDt {
    Float_t slope_;
    Float_t offset_;
    Float_t slopeUnc_;
    Float_t chisq_;
    Int_t   dof_;
    Float_t slope    () const { return slope_    ; }
    Float_t offset   () const { return offset_   ; }
    Float_t slopeUnc () const { return slopeUnc_ ; }
    Float_t chisq    () const { return chisq_    ; }
    Int_t   dof      () const { return dof_      ; }
    KalSeedDtDt() : slope_(0.f), offset_(0.f), slopeUnc_(0.f), chisq_(0.f), dof_(0) {}
  };

  using KalSeedDtDtCollection = std::vector<mu2e::KalSeedDtDt>;
}
#endif
