#ifndef CrvHitInfoMC_hh
#define CrvHitInfoMC_hh

#include "CLHEP/Vector/ThreeVector.h"
#include <vector>
#include "Rtypes.h"

namespace mu2e
{
  struct CrvHitInfoMC  //information about the MC track which most likely caused the CRV coincidence triplets
  {
    Bool_t              _valid;            //was an MC particle found that matches the coincidence triplets? 
    Int_t               _pdgId;            //PDG ID of this MC particle
    Int_t               _primaryPdgId;     //PDG ID of the primary particle of this MC particle (helps to determine whether it was a cosmic ray, etc.)
    Float_t             _primaryE;         //energy of the primary particle of this MC particle
    Float_t             _primaryX, _primaryY, _primaryZ;   //starting point of the primary particle of this MC particle
    Float_t             _x, _y, _z;        //position of the MC particle when it "created" the first StepPointMC
    Float_t             _time;             //time of the MC particle when it "created" the first StepPointMC
    Float_t             _depositedEnergy;  //total energy deposited for this cluster (not just for this track)
    CrvHitInfoMC(bool valid, int pdgId, 
              int primaryPdgId, float primaryE, CLHEP::Hep3Vector primaryPos, 
              CLHEP::Hep3Vector pos, float time, float depositedEnergy) :
              _valid(valid),
              _pdgId(pdgId),
              _primaryPdgId(primaryPdgId),
              _primaryE(primaryE),
              _primaryX(primaryPos.x()),
              _primaryY(primaryPos.y()),
              _primaryZ(primaryPos.z()),
              _x(pos.x()), _y(pos.y()), _z(pos.z()),
              _time(time),
              _depositedEnergy(depositedEnergy)
              {}
    CrvHitInfoMC() :
              _valid(false),
              _pdgId(0),
              _primaryPdgId(0),
              _primaryE(0),
              _primaryX(0),
              _primaryY(0),
              _primaryZ(0),
              _x(0), _y(0), _z(0),
              _time(0),
              _depositedEnergy(0)
              {}
  };

  typedef std::vector<CrvHitInfoMC>   CrvHitInfoMCCollection;    //this is the MC vector which will be stored in the main TTree 

}
#endif


