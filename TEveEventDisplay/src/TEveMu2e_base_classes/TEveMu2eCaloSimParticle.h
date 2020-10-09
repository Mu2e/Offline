#ifndef TEveMu2eCaloSimParticle_h
#define TEveMu2eCaloSimParticle_h

#include <TObject.h>
#include <string.h>

#include <TEvePointSet.h>
#include <TEveLine.h>
#include "MCDataProducts/inc/CaloHitSimPartMC.hh"

namespace mu2e {
  class TEveMu2eCaloSimParticle : public TEvePointSet {
    public:
      #ifndef __CINT__
      explicit TEveMu2eCaloSimParticle();
      TEveMu2eCaloSimParticle(CaloHitSimPartMC simpartmc) : fCaloHitSimPartMC(simpartmc){};
      virtual ~TEveMu2eCaloSimParticle(){};

      CaloHitSimPartMC fCaloHitSimPartMC; 
      Int_t mSize= 1; 


      void DrawParticle3D(const std::string &pstr, Int_t b,CLHEP::Hep3Vector pointInMu2e, TEveElementList *ParticleList);
      void DrawParticle2D(const std::string &pstr, Int_t b,CLHEP::Hep3Vector pointInMu2e, TEveElementList *ParticleList);



      #endif
      ClassDef(TEveMu2eCaloSimParticle, 0);
    };
    typedef std::vector<mu2e::TEveMu2eCaloSimParticle> TEveMu2eCaloSimParticleCollection;
}
#endif
