#ifndef TrkDiag_CosmicTrackInfo_HH
#define TrkDiag_CosmicTrackInfo_HH

#include "RecoDataProducts/inc/XYZVec.hh"
#include "Rtypes.h"
namespace mu2e {
  struct CosmicTrackInfo {
    Int_t N_hits;
    Int_t ChiSq_NDF;
    Int_t ChiSq;
    Double_t M_XY;
    Double_t C_XY;
    
    
  };

 
#endif
