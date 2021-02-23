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
      TEveMu2eCustomHelix();
      explicit TEveMu2eCustomHelix(const TEveMu2eCustomHelix &helix);
      explicit TEveMu2eCustomHelix(HelixSeed  const& hseed);
      explicit TEveMu2eCustomHelix(KalSeed kseed);
      virtual ~TEveMu2eCustomHelix(){};
      #endif
      
      KalSeed fKalSeed; 
      HelixSeed fHelixSeed;
      TrkExtTraj fTrkExtTraj;

      void DrawHelixTrack();
      void Draw2DProjection();
      //TODO - this class will need redesigning
      void SetSeedInfo(KalSeed const&  seed);
      void SetPostionAndDirectionFromHelixSeed(double zpos);
      void SetPostionAndDirectionFromKalRep(double zpos);
      void SetMomentumExt();
      void SetParticleExt();
      const std::string Title();
      std::string DataTitle(const std::string &pstr, int n);
      XYZVec Direction;
      XYZVec Position;
      double Momentum;
      int PDGcode;
      double Charge;
      double Mass;
      bool _trajectory;
      ClassDef( TEveMu2eCustomHelix, 0);
    };
}
#endif
