#ifndef CosmicRayShieldGeom_CRSScintillatorLayer_hh
#define CosmicRayShieldGeom_CRSScintillatorLayer_hh
//
// Representation of one Scintillator Layer in  CosmicRayShield
//
//
// Original author KLG; somewhat based on  Rob Kutschke's Layer
//

#include <vector>

#include "Offline/CosmicRayShieldGeom/inc/CRSScintillatorLayerId.hh"
#include "Offline/CosmicRayShieldGeom/inc/CRSScintillatorBar.hh"

#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e
{

  class CRSScintillatorLayer
  {

    friend class CRSScintillatorModule;
    friend class CRSScintillatorShield;
    friend class CosmicRayShieldMaker;

    CRSScintillatorLayer();

    public:

    CRSScintillatorLayer(CRSScintillatorLayerId const & id);

    // Accept the compiler generated destructor, copy constructor and assignment operators

    CRSScintillatorLayerId const & id() const { return _id;}

    int nBars() const { return int(_bars.size()); }

    CRSScintillatorBar const & getBar( int n ) const
    {
      return *_bars.at(n);
    }

    CRSScintillatorBar const & getBar( const CRSScintillatorBarId& id ) const
    {
      return getBar(id.getBarNumber());
    }

    const std::vector<std::shared_ptr<CRSScintillatorBar> >& getBars() const { return _bars; }

    const CLHEP::Hep3Vector &getPosition() const
    {return _position;}

    const std::vector<double> &getHalfLengths() const {return _halfLengths;}
    double getHalfThickness() const { return _halfLengths[_localToWorld[0]];}
    double getHalfWidth() const { return _halfLengths[_localToWorld[1]];}
    double getHalfLength() const { return _halfLengths[_localToWorld[2]];}

    // Formatted string embedding the id of the layer.
    std::string name( std::string const & base ) const;

    private:

    CRSScintillatorLayerId _id;

    // Pointers to the bars in this layer.
    // They refer to the same objects as the bars in CosmicRayShield (which holds all CRV bars)
    std::vector<std::shared_ptr<CRSScintillatorBar> > _bars;

    CLHEP::Hep3Vector _position;
    std::vector<double> _halfLengths;
    std::vector<int> _localToWorld;  //0th entry: thickness
                                     //1st entry: width
                                     //2nd entry: length
  };
}

#endif /* CosmicRayShieldGeom_CRSScintillatorLayer_hh */
