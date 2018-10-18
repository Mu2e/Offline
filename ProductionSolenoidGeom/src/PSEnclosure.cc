#include "ProductionSolenoidGeom/inc/PSEnclosure.hh"

namespace mu2e {

  // genreflex perstency requires default ctr
  PSEnclosure::PSEnclosure() {}

  std::ostream& operator<<(std::ostream& os, const PSEnclosure& pse) {
    if ( pse.version() == 1 ) {
      os<<"PSEnclosure("
	<<"material="<<pse.shell().materialName()
	<<", OD="<<2*pse.shell().outerRadius()
	<<", length="<<2*(pse.shell().halfLength() + pse.endPlate().halfLength())
	<<", endPlate.thickness="<<2*pse.endPlate().halfLength()
	<<", shell.thickness="<<(pse.shell().outerRadius() - pse.shell().innerRadius())
	<<", windows = { "
	;

      for(unsigned i=0; i<pse.nWindows(); ++i) {
	const CLHEP::Hep3Vector tmp(pse.windows()[i].originInMu2e() - pse.endPlate().originInMu2e());
	os<<"Window("
	  <<"materialName="<<pse.windows()[i].materialName()
	  <<", thickness="<<2*pse.windows()[i].halfLength()
	  <<", r="<<pse.windows()[i].outerRadius()
	  <<", xoff="<<tmp.x()
	  <<", yoff="<<tmp.y()
	  <<"), ";

      }

    } else {

      os<<"PSEnclosure("
	<<"material="<<pse.shellCone().materialName()
	<<", OD1="<<2*pse.shellCone().outerRadius1()
	<<", OD2="<<2*pse.shellCone().outerRadius2()
	<<", length="<<2*(pse.shellCone().halfLength() + pse.endPlate().halfLength())
	<<", endPlate.thickness="<<2*pse.endPlate().halfLength()
	<<", shell.thickness="<<(pse.shellCone().outerRadius1() - pse.shellCone().innerRadius1())
	<<", windows = { "
	;

      for(unsigned i=0; i<pse.nWindows(); ++i) {
	const CLHEP::Hep3Vector tmp(pse.windows()[i].originInMu2e() - pse.endPlate().originInMu2e());
	os<<"Window("
	  <<"materialName="<<pse.windows()[i].materialName()
	  <<", thickness="<<2*pse.windows()[i].halfLength()
	  <<", r="<<pse.windows()[i].outerRadius()
	  <<", xoff="<<tmp.x()
	  <<", yoff="<<tmp.y()
	  <<"), ";

      }

    }

    os<<" } )";
    return os;
  }

} // namespace mu2e
