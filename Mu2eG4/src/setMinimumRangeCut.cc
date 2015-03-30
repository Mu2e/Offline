//
// Set the G4 minimum range cut as specified in the geometry file.
//
// $Id: setMinimumRangeCut.cc,v 1.2 2012/07/15 22:06:17 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/07/15 22:06:17 $
//

#include "Mu2eG4/inc/setMinimumRangeCut.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "fhiclcpp/ParameterSet.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "G4UImanager.hh"

#include <string>
#include <sstream>

namespace mu2e{

  void setMinimumRangeCut(double minRangeCut){
    std::ostringstream out;
    out << "/run/setCut " << minRangeCut << " mm ";
    G4UImanager* UI = G4UImanager::GetUIpointer();
    UI->ApplyCommand(out.str().c_str());
    mf::LogInfo("GEOM")
      << "Setting minRange cut to " << minRangeCut << " mm\n";
  }

  void setMinimumRangeCut( SimpleConfig const& config ){
    // If this parameter is absent, leave the default from G4 unchanged.
    std::string name("g4.minRangeCut");
    if ( config.hasName(name) ){
      double minRangeCut = config.getDouble(name);
      setMinimumRangeCut(minRangeCut);
    }
  }

  void setMinimumRangeCut(const fhicl::ParameterSet& pset){
    setMinimumRangeCut(pset.get<double>("physics.minRangeCut"));
  }

}  // end namespace mu2e
