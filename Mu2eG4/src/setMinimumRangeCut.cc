//
// Set the G4 minimum range cut as specified in the geometry file.
//
// $Id: setMinimumRangeCut.cc,v 1.1 2012/06/04 19:28:01 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/06/04 19:28:01 $
//

#include "Mu2eG4/inc/setMinimumRangeCut.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "G4UImanager.hh"

#include <string>
#include <sstream>

namespace mu2e{

  void setMinimumRangeCut( SimpleConfig const& config ){

    // If this parameter is absent, leave the default from G4 unchanged.
    std::string name("g4.minRangeCut");
    if ( config.hasName(name) ){
      double minRangeCut = config.getDouble(name);
      std::ostringstream out;
      out << "/run/setCut " << minRangeCut << " mm ";
      G4UImanager* UI = G4UImanager::GetUIpointer();
      UI->ApplyCommand(out.str().c_str());
      mf::LogInfo("GEOM")
        << "Setting minRange cut to " << minRangeCut << " mm\n";
    }
  }

}  // end namespace mu2e
