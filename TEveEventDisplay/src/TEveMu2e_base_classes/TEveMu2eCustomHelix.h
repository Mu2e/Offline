#ifndef TEveMu2eCustomHelix_h
#define TEveMu2eCustomHelix_h

#include <TObject.h>
#include <THelix.h>
#include <TEveLine.h>
#include "RecoDataProducts/inc/HelixSeed.hh"
#include "RecoDataProducts/inc/KalSeed.hh"
#include "RecoDataProducts/inc/TrkExtTraj.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"

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
      TrkExtTraj fTrkExtTraj_;

      void DrawHelixTrack();
      void Draw2DProjection();

      void SetSeedInfo(KalSeed seed);
      void SetPostionAndDirectionFromHelixSeed(double zpos);
      void SetPostionAndDirectionFromKalRep(double zpos);
      void SetMomentumExt();
      void SetParticleExt();

      XYZVec Direction_;
      XYZVec Position_;
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
