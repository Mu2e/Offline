#ifndef CosmicRayShieldGeom_CRSScintillatorBarDetail_hh
#define CosmicRayShieldGeom_CRSScintillatorBarDetail_hh

//
// Representation of common properties of the Scintillator Bars etc...
//
// $Id: CRSScintillatorBarDetail.hh,v 1.6 2013/09/13 06:42:44 ehrlich Exp $
// $Author: ehrlich $
// $Date: 2013/09/13 06:42:44 $
//
// Original author KLG; somewhat based on Rob Kutschke's StrawDetail
//

#include <string>
#include <vector>

namespace mu2e 
{

  class CRSScintillatorBarDetail
  {

    friend class CosmicRayShieldMaker;

    public:

    CRSScintillatorBarDetail(): _materialName(), _halfLengths() {}

    CRSScintillatorBarDetail(std::string const & materialName,
                             std::vector<double> const & halfLengths);

    // Compiler generated versions are OK for destructor
    // and for copy and assignment constructors.

    std::string const & getMaterialName() const { return _materialName;}

    std::vector<double> const & getHalfLengths() const { return _halfLengths;}

    private:

    std::string _materialName;
    std::vector<double> _halfLengths;
  };

}  //namespace mu2e

#endif /* CosmicRayShieldGeom_CRSScintillatorBarDetail_hh */
