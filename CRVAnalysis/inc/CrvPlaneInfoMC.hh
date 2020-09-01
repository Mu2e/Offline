#ifndef CrvPlaneInfoMC_hh
#define CrvPlaneInfoMC_hh

#include "CLHEP/Vector/ThreeVector.h"
#include <vector>
#include "Rtypes.h"

namespace mu2e
{
  struct CrvPlaneInfoMC  //information about the point where the MC trajectory crosses the xz plane of CRV-T
  {
    Int_t               _pdgId;            //PDG ID of this MC particle
    Int_t               _primaryPdgId;     //PDG ID of the primary particle of this MC particle
    Float_t             _primaryE;         //energy of the primary particle of this MC particle
    Float_t             _primaryX, _primaryY, _primaryZ;   //starting point of the primary particle of this MC particle
    Float_t             _x, _y, _z;        //position of the MC particle when it crosses the xz plane of CRV-T
    Float_t             _xDir, _yDir, _zDir;    //direction of the MC particle when it crosses the xz plane of CRV-T
    Float_t             _time;             //time of the MC particle when it crosses the xz plane of CRV-T
    Float_t             _kineticEnergy;    //time of the MC particle when it crosses the xz plane of CRV-T
    Int_t               _dataSource;       //temporary variable; will be removed (1...data from stepPointMCs, 2...data from trajectory extrapolation)
    CrvPlaneInfoMC(int pdgId, int primaryPdgId, float primaryE, CLHEP::Hep3Vector primaryPos, 
              CLHEP::Hep3Vector pos, CLHEP::Hep3Vector dir, float time, float kineticEnergy, int dataSource) :
              _pdgId(pdgId),
              _primaryPdgId(primaryPdgId),
              _primaryE(primaryE),
              _primaryX(primaryPos.x()),
              _primaryY(primaryPos.y()),
              _primaryZ(primaryPos.z()),
              _x(pos.x()), _y(pos.y()), _z(pos.z()),
              _xDir(dir.x()), _yDir(dir.y()), _zDir(dir.z()),
              _time(time),
              _kineticEnergy(kineticEnergy),
              _dataSource(dataSource)
              {}
    CrvPlaneInfoMC() :
              _pdgId(0),
              _primaryPdgId(0),
              _primaryE(0),
              _primaryX(0),
              _primaryY(0),
              _primaryZ(0),
              _x(0), _y(0), _z(0),
              _xDir(0), _yDir(0), _zDir(0),
              _time(0),
              _kineticEnergy(0),
              _dataSource(0)
              {}
  };

  typedef std::vector<CrvPlaneInfoMC>   CrvPlaneInfoMCCollection;    //this is the MC vector which will be stored in the main TTree 

}
#endif


