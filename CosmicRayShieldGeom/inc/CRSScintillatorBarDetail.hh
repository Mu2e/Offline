#ifndef CosmicRayShieldGeom_CRSScintillatorBarDetail_hh
#define CosmicRayShieldGeom_CRSScintillatorBarDetail_hh

//
// Representation of common properties of the Scintillator Bars etc...
//
// $Id: CRSScintillatorBarDetail.hh,v 1.5 2011/05/20 20:21:47 wb Exp $
// $Author: wb $
// $Date: 2011/05/20 20:21:47 $
//
// Original author KLG; somewhat based on Rob Kutschke's StrawDetail
//

#include <string>
#include <vector>

namespace mu2e {

  class CRSScintillatorBarDetail{

    friend class CosmicRayShieldMaker;

  public:

    CRSScintillatorBarDetail():
      _id(-1),
      _materialNames(),
      _halfLengths()
    {}

    CRSScintillatorBarDetail( int const id,
                              std::vector<std::string> const & materialNames,
                              std::vector<double> const & halfLengths
                              );

    // Compiler generated versions are OK for destructor
    // and for copy and assignment constructors.

    int const  Id() const { return  _id; }

    std::string const & getMaterialName(int idx) const { return _materialNames.at(idx);}

    std::vector<std::string> const & getMaterialNames() const { return _materialNames;}

    std::vector<double> const & getHalfLengths() const { return _halfLengths;}

    std::string name( std::string const & base ) const;

  private:

    // Identifier for this type of CRSScintillatorBar.
    int _id;

    // Order of materials is:

    // _scintillatorBarMaterialName
    // _scintillatorBarPigmentationMaterialName
    // _scintillatorModuleOuterSheetMaterialName
    // _scintillatorModuleInterLayerSheetMaterialName

    std::vector<std::string> _materialNames;

    // outer dimensions
    std::vector<double> _halfLengths;

  };

}  //namespace mu2e

#endif /* CosmicRayShieldGeom_CRSScintillatorBarDetail_hh */
