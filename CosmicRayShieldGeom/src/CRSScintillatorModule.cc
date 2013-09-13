//
// Representation of one Scintillator Module in  CosmicRayShield
//
//
// $Id: CRSScintillatorModule.cc,v 1.5 2013/09/13 06:42:44 ehrlich Exp $
// $Author: ehrlich $
// $Date: 2013/09/13 06:42:44 $
//
// Original author KLG based on Rob Kutschke's Sector
//

#include <sstream>

#include "CosmicRayShieldGeom/inc/CRSScintillatorModule.hh"

using namespace std;

namespace mu2e 
{

    CRSScintillatorModule::CRSScintillatorModule():
      _id(CRSScintillatorModuleId(-1,-1))
    {}

    CRSScintillatorModule::CRSScintillatorModule( CRSScintillatorModuleId const& id):
      _id(id)
    {}

    string CRSScintillatorModule::name( string const& base ) const
    {
      ostringstream os;

      os << base
         << _id.getShieldNumber() << "_"
         << _id.getModuleNumber();
      return os.str();
    }

} // end namespace mu2e
