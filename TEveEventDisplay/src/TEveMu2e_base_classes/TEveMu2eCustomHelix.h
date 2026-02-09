#ifndef TEveMu2eCustomHelix_h
#define TEveMu2eCustomHelix_h

#include <TObject.h>
#include <THelix.h>
#include <TEveLine.h>
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"

using namespace mu2e;

namespace mu2e {
  class   TEveMu2eCustomHelix: public TEveLine {
    public:
      #ifndef __CINT__
      explicit  TEveMu2eCustomHelix();
      TEveMu2eCustomHelix(const TEveMu2eCustomHelix &helix);
      virtual ~TEveMu2eCustomHelix(){};
      #endif

      KalSeed fKalSeed_;
      HelixSeed fHelixSeed_;

      void DrawHelixTrack();
      void Draw2DProjection();

      void SetSeedInfo(KalSeed seed);

      XYZVectorF Direction_;
      XYZVectorF Position_;
      double Momentum_;
      int PDGcode_;
      double Charge_;
      double Mass_;
      double Time_;
      double Radius_;
      ClassDef( TEveMu2eCustomHelix, 0);
    };
}
#endif
