#include "GeometryService/inc/ProductionTargetMaker.hh"

#include "cetlib/exception.h"

#include "CLHEP/Units/SystemOfUnits.h"

#include "Mu2eUtilities/inc/SimpleConfig.hh"

#include "GeometryService/inc/ProductionTarget.hh"

namespace mu2e {

  ProductionTargetMaker::ProductionTargetMaker(const SimpleConfig& c) 
    : m_det(new ProductionTarget(c.getDouble("targetPS_rOut"),
				 c.getDouble("targetPS_halfLength"),
				 c.getDouble("targetPS_rotX") * CLHEP::degree,
				 c.getDouble("targetPS_rotY") * CLHEP::degree,
				 c.getHep3Vector("productionTarget.position")
				 )
	    )
  {}

}
