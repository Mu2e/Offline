//
// Representation of one CRSScintillatorBar in CosmicRayShield
//
// $Id: CRSScintillatorBar.cc,v 1.4 2013/09/13 06:42:44 ehrlich Exp $
// $Author: ehrlich $
// $Date: 2013/09/13 06:42:44 $
//
// Original author KLG
//

#include <sstream>
#include "CosmicRayShieldGeom/inc/CRSScintillatorBar.hh"
#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"

namespace mu2e 
{
    CRSScintillatorBar::CRSScintillatorBar(CRSScintillatorBarIndex const &index, 
                       CRSScintillatorBarId const &id,
                       CLHEP::Hep3Vector const &position,
                       const std::shared_ptr<CRSScintillatorBarDetail> detail) : 
    _index(index),
    _id(id),
    _position(position),
    _detail(detail)
    {
    }

    std::string CRSScintillatorBar::name( std::string const& base ) const
    {
      std::ostringstream os;
      os << base
       << _id.getShieldNumber() << "_"
       << _id.getModuleNumber() << "_"
       << _id.getLayerNumber()  << "_"
       << _id.getBarNumber();
      return os.str();
    }

} // namespace mu2e

