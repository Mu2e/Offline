#ifndef TrackerConditions_StrawAlignment_hh
#define TrackerConditions_StrawAlignment_hh
//
// StrawAlignment describes individual straws as extended objects in space
// Original author: Dave Brown (LBNL) 7/22
//
#include "TrackerGeom/inc/Straw.hh"
#include "TrackerConditions/inc/StrawEndAlignment.hh"
#include "TrackerConditions/inc/StrawTension.hh"
#include "Mu2eInterfaces/inc/ProditionsEntity.hh"
#include "DataProducts/inc/StrawId.hh"
#include "Math/Vector3D.h"
#include <iostream>
#include <map>

namespace mu2e {

  class StrawAlignment : virtual public ProditionsEntity {
    public:
      typedef ROOT::Math::XYZVector Vec3; // spatial only vector
      typedef std::shared_ptr<StrawAlignment> ptr_t;
      typedef std::shared_ptr<const StrawAlignment> cptr_t;

      explicit StrawAlignment() : name_("StrawAlignment") {}

      // describe the position (in local UVW coordinates) of the straw and wire centers at a specified length
      // This uses the proditions alignment info, plus the High Voltage and the direction of gravity
      Vec3 wireCenterUVW(Straw const& straw, double Upos, double HV, Vec3 gravity=Vec3(0.0,0.0,9800)) const;
      Vec3 strawCenterUVW(Straw const& straw, double Upos, double HV, Vec3 gravity=Vec3(0.0,0.0,9800)) const;

      // XYZ position of the straw and wire centers at a specified length
      Vec3 wireCenter(Straw const& straw, double Upos, double HV, Vec3 gravity=Vec3(0.0,0.0,9800)) const;
      Vec3 strawCenter(Straw const& straw, double Upos, double HV, Vec3 gravity=Vec3(0.0,0.0,9800)) const;


  private:
    Vec3 centerUVW(StrawEndAlignment const& sea, Straw const& straw, double Upos, double HV, Vec3 gravity=Vec3(0.0,0.0,9800)) const;
    std::string name_; // conditions name
    std::map<StrawId, StrawEndAlignment> wireEnds_; // describe the wire ends
    std::map<StrawId, StrawEndAlignment> strawEnds_; // describe the (mylar shell) straw ends
    std::map<StrawId, StrawTension> wireTension_; // wire tensions
    std::map<StrawId, StrawTension> strawTension_; // straw tensions
  };
}
#endif
