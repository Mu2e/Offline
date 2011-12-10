//
// Make a ExtinctionMonitor.
//
// $Id: ExtMonUCI.hh,v 1.1 2011/12/10 00:16:15 youzy Exp $
// $Author: youzy $
// $Date: 2011/12/10 00:16:15 $

#ifndef EXTMONUCI_HH
#define EXTMONUCI_HH

#include <vector>

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"

#include "GeometryService/inc/Detector.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/Mu2eBuilding.hh"

#include "ExtinctionMonitorUCIGeom/inc/ExtMonUCICol.hh"

using namespace std;

namespace mu2e {

  namespace ExtMonUCI {
    
    class ExtMon : public Detector {

    friend class ExtMonMaker;

    public: 
      // implement Detector's method
      virtual std::string name() const { return "ExtMonUCI"; }

      // envelope
      const vector<double>& envelopeParams() const { return _envelopeParams; }
      const vector<double>& envelopeOrigin() const { return _envelopeOrigin; }

      void setHallOriginInMu2e(double x, double y, double z) { _hallOriginInMu2e = CLHEP::Hep3Vector(x,y,z); }

      CLHEP::Hep3Vector origin() const;
      CLHEP::Hep3Vector originLocal() const;
      const CLHEP::HepRotation& rotation() const { return _rotation; }
      const CLHEP::Hep3Vector&  hallOriginInMu2e() const { return _hallOriginInMu2e; }

      // collimator
      int nCols() const { return _nCols; }
      const vector<double>& colOuterHalfLengths() const { return _colOuterHalfLengths; }
      const vector<double>& colInnerHalfLengths() const { return _colInnerHalfLengths; }
      const vector<double>& colPosition1() const { return _colPosition1; }
      const vector<double>& colPosition2() const { return _colPosition2; }
      const ExtMonCol* col(unsigned int iCol) const; 

      // Coordinate conversion to/from the Mu2e frame
      // 
      // - The (0,0,0) point is in the middle of the ExtMon detector.
      // 
      // In G4 this is the coordinate system of the ExtMonUCI volume.
      CLHEP::Hep3Vector extMonToMu2ePoint( CLHEP::Hep3Vector const& v ) const;      
      CLHEP::Hep3Vector extMonToMu2eMomentum( CLHEP::Hep3Vector const& v ) const;      
      CLHEP::Hep3Vector mu2eToExtMonPoint( CLHEP::Hep3Vector const& v ) const;
      CLHEP::Hep3Vector mu2eToExtMonMomentum( CLHEP::Hep3Vector const& v ) const;

    private: 

      ExtMon(const std::vector<double>& envelopeParams, 
	     const std::vector<double>& envelopeOrigin);

      CLHEP::Hep3Vector   _hallOriginInMu2e;

      std::vector<double> _envelopeParams;
      std::vector<double> _envelopeOrigin;
      CLHEP::HepRotation  _rotation;

      std::vector<double> _abs1Params;
      std::vector<double> _abs1Origin;

      // Collimators
      int _nCols;
      std::vector<double> _colOuterHalfLengths;
      std::vector<double> _colInnerHalfLengths;
      std::vector<double> _colPosition1;
      std::vector<double> _colPosition2;
      std::vector<ExtMonCol> _cols;
    };

  }
}

#endif/*EXTMONUCI_HH*/
