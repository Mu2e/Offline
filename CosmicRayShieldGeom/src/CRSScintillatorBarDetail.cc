//
// Representation of common properties of the Scintillator Bars etc...
//
// $Id: CRSScintillatorBarDetail.cc,v 1.5 2013/09/13 06:42:44 ehrlich Exp $
// $Author: ehrlich $
// $Date: 2013/09/13 06:42:44 $
//
// Original author KLG; somewhat based on Rob Kutschke StrawDetail
//

#include <sstream>

#include "CosmicRayShieldGeom/inc/CRSScintillatorBarDetail.hh"

using namespace std;

namespace mu2e 
{
  CRSScintillatorBarDetail::CRSScintillatorBarDetail(std::string const& materialName,
                                                     std::vector<double> const& halfLengths) :
    _materialName(materialName),
    _halfLengths(halfLengths)
  {}

} // namespace mu2e
