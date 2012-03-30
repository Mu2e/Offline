#ifndef ProductionSolenoidGeom_PSEnclosureMaker_hh
#define ProductionSolenoidGeom_PSEnclosureMaker_hh
//
// Class to construct and return PSEnclosure
//
// $Id: PSEnclosureMaker.hh,v 1.2 2012/03/30 04:14:38 gandr Exp $
// $Author: gandr $
// $Date: 2012/03/30 04:14:38 $
//
// Original author Andrei Gaponenko
//

#include <memory>
#include <string>
#include "boost/utility.hpp"

namespace CLHEP { class Hep3Vector; }

namespace mu2e {

  class PSEnclosure;
  class SimpleConfig;

  class PSEnclosureMaker : boost::noncopyable {

  public:

    PSEnclosureMaker(const SimpleConfig& config,

                     // The center of the downstream surface of the PS
                     const CLHEP::Hep3Vector& psEndRefPoint,

                     const std::string& psInsideMaterialName);

    std::auto_ptr<PSEnclosure> getPtr() { return pse_; }

  private:

    std::auto_ptr<PSEnclosure> pse_;

  };

}  //namespace mu2e

#endif /* ProductionSolenoidGeom_PSEnclosureMaker_hh */
