#ifndef TEveMu2eCustomHelix_h
#define TEveMu2eCustomHelix_h

#include <TObject.h>
#include <THelix.h>
#include <TEveLine.h>
#include "RecoDataProducts/inc/HelixSeed.hh"
#include "RecoDataProducts/inc/KalSeed.hh"
#include "RecoDataProducts/inc/TrkExtTraj.hh"

using namespace mu2e;

namespace mu2e {
  class   TEveMu2eCustomHelix: public TEveLine {
    public:
      #ifndef __CINT__
      explicit  TEveMu2eCustomHelix();
      TEveMu2eCustomHelix(const TEveMu2eCustomHelix &helix);
      TEveMu2eCustomHelix(HelixSeed hseed);
      TEveMu2eCustomHelix(KalSeed kseed);
      virtual ~TEveMu2eCustomHelix(){};
      #endif
      
      KalSeed fKalSeed; 
      HelixSeed fHelixSeed;
      TrkExtTraj fTrkExtTraj;

      void DrawHelixTrack();
      void Draw2DProjection();

      void SetSeedInfo(KalSeed seed);
      void SetPostionAndDirectionFromHelixSeed(double zpos);
      void SetPostionAndDirectionFromKalRep(double zpos);
      void SetMomentumExt();
      void SetParticleExt();

      XYZVec Direction;
      XYZVec Position;
      double Momentum;
      int PDGcode;
      double Charge;
      double Mass;
      double Time;
      double Radius;
      bool _trajectory;
      ClassDef( TEveMu2eCustomHelix, 0);
    };
}
#endif
