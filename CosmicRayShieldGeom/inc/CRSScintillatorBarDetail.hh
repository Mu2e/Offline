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
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e 
{

  class CRSScintillatorBarDetail
  {

//    friend class CosmicRayShieldMaker;

    //disable default constructor
    CRSScintillatorBarDetail();

    public:

    CRSScintillatorBarDetail(std::string const & materialName,
                             std::vector<double> const & halfLengths,
                             std::vector<int> const & localToWorld,
                             std::string const & CMBmaterialName,
                             double CMBoffset, double CMBhalfThickness,
                             bool CMBside0, bool CMBside1,
                             double fiberSeparation);

    // Compiler generated versions are OK for destructor
    // and for copy and assignment constructors.

    std::string const & getMaterialName() const { return _materialName;}

    std::vector<double> const & getHalfLengths() const { return _halfLengths;}

    double getHalfThickness() const { return _halfLengths[_localToWorld[0]];}
    double getHalfWidth() const { return _halfLengths[_localToWorld[1]];}
    double getHalfLength() const { return _halfLengths[_localToWorld[2]];}

    CLHEP::Hep3Vector toWorld(const CLHEP::Hep3Vector &localPosition, const CLHEP::Hep3Vector &barPosition) const;
    CLHEP::Hep3Vector toLocal(const CLHEP::Hep3Vector &worldPosition, const CLHEP::Hep3Vector &barPosition) const;
    CLHEP::Hep3Vector toLocalNormalized(const CLHEP::Hep3Vector &worldPosition, const CLHEP::Hep3Vector &barPosition) const;
    bool isInside(const CLHEP::Hep3Vector &worldPosition, const CLHEP::Hep3Vector &barPosition) const;

    int getThicknessDirection() const {return _localToWorld[0];}
    int getWidthDirection() const {return _localToWorld[1];}
    int getLengthDirection() const {return _localToWorld[2];}

    private:

    std::string _materialName;
    std::vector<double> _halfLengths;
    std::vector<int> _localToWorld;  //0th entry: thickness
                                     //1st entry: width
                                     //2nd entry: length


    /********************/
    // counter motherboard section

    public:

    std::string const & getCMBMaterialName() const { return _CMBmaterialName;}

    double getCMBOffset() const {return _CMBoffset;}
    double getCMBHalfThickness() const {return _CMBhalfThickness;}

    CLHEP::Hep3Vector getCMBPosition(int side, const CLHEP::Hep3Vector &barPosition) const;  //"side" can be 0 or 1

    std::vector<double> getCMBHalfLengths() const;

    bool hasCMB(int side) const;

    private:

    std::string _CMBmaterialName;
    double _CMBoffset;
    double _CMBhalfThickness;
    bool   _CMBside0;
    bool   _CMBside1;

    /********************/
    // SiPM section

    public:

    double getFiberSeparation() const {return _fiberSeparation;}

    CLHEP::Hep3Vector getSiPMPosition(int SiPMNumber, const CLHEP::Hep3Vector &barPosition) const;
    CLHEP::Hep3Vector getSiPMPosition(int fiberNumber, int side, const CLHEP::Hep3Vector &barPosition) const;

    private:
 
    double _fiberSeparation;

  };

}  //namespace mu2e

#endif /* CosmicRayShieldGeom_CRSScintillatorBarDetail_hh */
