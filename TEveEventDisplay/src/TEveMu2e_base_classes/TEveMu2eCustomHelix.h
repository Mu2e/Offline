#ifndef TEveMu2eCustomHelix_h
#define TEveMu2eCustomHelix_h

#include <TObject.h>
#include <THelix.h>
#include <TPolyLine3D.h>
#include "RecoDataProducts/inc/HelixSeed.hh"
#include "RecoDataProducts/inc/KalSeed.hh"

#include <TEveTrackPropagator.h>

namespace mu2e {
  class   TEveMu2eCustomHelix: public TEveLine {
    public:
       #ifndef __CINT__
      explicit  TEveMu2eCustomHelix(){};
      TEveMu2eCustomHelix(const TEveMu2eCustomHelix &helix) { fKalSeed = helix.fKalSeed;} ;
      TEveMu2eCustomHelix(HelixSeed hseed){fHelixSeed = hseed;};
      TEveMu2eCustomHelix(KalSeed kseed){
        fKalSeed = kseed;
        this->Momentum = fKalSeed.helix()->helix().momentum();
        this->PDGcode = fKalSeed.particle().particleType();
        this->Charge = fKalSeed.particle().charge();
        this->Mass = fKalSeed.particle().mass();
      };
      virtual ~ TEveMu2eCustomHelix(){};
      #endif
      
      KalSeed fKalSeed; 
      HelixSeed fHelixSeed;
      TrkExtTraj fTrkExtTraj;

      void DrawHelixTrack();
      void Draw2DProjection();

      void SetSeedInfo(KalSeed seed) { 
        fKalSeed = seed;
        this->Momentum = fKalSeed.helix()->helix().momentum();
        this->PDGcode = fKalSeed.particle().particleType();
        this->Charge = fKalSeed.particle().charge();
        this->Mass = fKalSeed.particle().mass();
      }

      void SetPostionAndDirectionFromHelixSeed(double zpos){
        fHelixSeed.helix().position(Position);
        fHelixSeed.helix().direction(zpos, Direction);
      }

      void SetPostionAndDirectionFromKalRep(double zpos){
        fKalSeed.helix()->helix().position(Position);
        fKalSeed.helix()->helix().direction(zpos, Direction);
      }

      void SetMomentumExt(){
        this->Momentum = fTrkExtTraj.front().momentum().mag();
      }

      void SetParticleExt(){
        this->PDGcode = 11; //FIXME
      }

      const std::string Title(){
        const std::string title = "Track PDG " + to_string(PDGcode) +" Momentum = " + to_string(Momentum) + " Charge "+ to_string(Charge);
        return title;
      }

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
