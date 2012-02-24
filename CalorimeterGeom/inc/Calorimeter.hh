#ifndef CalorimeterGeom_Calorimeter_hh
#define CalorimeterGeom_Calorimeter_hh
//
// Hold all geometry and identifier information about
// a Calorimeter.  In order to insulate this class from
// knowledge of databases etc, this class must not know
// how to make itself.
//
// $Id: Calorimeter.hh,v 1.15 2012/02/24 20:55:48 gandr Exp $
// $Author: gandr $
// $Date: 2012/02/24 20:55:48 $
//
// Original author R. Bernstein and Rob Kutschke
//

#include <vector>

//
// Mu2e includes
#include "Mu2eInterfaces/inc/Detector.hh"
#include "CalorimeterGeom/inc/Vane.hh"

namespace mu2e {
    class Calorimeter: virtual public Detector{

      friend class CalorimeterMaker;

    public:
      Calorimeter(){}
      ~Calorimeter(){}

      CLHEP::Hep3Vector const& getOrigin() const { return _origin; }

      unsigned int nVane() const { return _nVane; }
      Vane const& getVane(int i) const { return _vanes.at(i); }

      unsigned int nCrystalPerVane() const { return _nCrystalZ*_nCrystalR; }
      unsigned int nCrystal()  const { return _nVane*_nCrystalZ*_nCrystalR; }
      unsigned int nCrystalR() const { return _nCrystalR; }
      unsigned int nCrystalZ() const { return _nCrystalZ; }

      unsigned int nROPerCrystal() const { return _nROPerCrystal; }
      unsigned int nRO() const { return _nROPerCrystal*_nVane*_nCrystalZ*_nCrystalR; }

      double crystalHalfSize()   const { return _crystalHW; }
      double crystalHalfLength() const { return _crystalHL; }
      double wrapperHalfThickness() const { return _wrapperHalfThickness; }
      double roHalfSize() const { return _roHalfTrans; }
      double roHalfThickness() const { return _roHalfThickness; }

      double getNonuniformity() const { return _nonUniformity; }
      double getTimeGap() const { return _timeGap; }
      double getElectronEdep() const { return _electronEdep; }
      double getElectronEmin() const { return _electronEmin; }

      // Vane ID (0..nVanes-1)
      int getVaneByRO(int roid) const {
        return roid/(_nCrystalZ*_nCrystalR*_nROPerCrystal);
      }
      // Crystal ID (0..Number_of_crystals_in_calorimeter-1)
      int getCrystalByRO(int roid) const { return (roid/_nROPerCrystal); }
      // Crystal ID within a vane (0..Number_of_crystals_in_vane-1)
      int getCrystalVaneByRO(int roid) const {
        return (roid/_nROPerCrystal)%(_nCrystalZ*_nCrystalR);
      }
      // Crystal R-coordinate within a vane (0..nCrystalR-1)
      int getCrystalRByRO(int roid) const {
        return ((roid/_nROPerCrystal)%(_nCrystalZ*_nCrystalR))/_nCrystalZ;
      }
      // Crystal Z-coordinate within a vane (0..nCrystalZ-1)
      int getCrystalZByRO(int roid) const {
        return ((roid/_nROPerCrystal)%(_nCrystalZ*_nCrystalR))%_nCrystalZ;
      }

      // Transfer Mu2e coordinates to local crystal coordinates
      CLHEP::Hep3Vector toCrystalFrame(int, CLHEP::Hep3Vector const&) const;

      // Get crystal center in Mu2e coordinates
      CLHEP::Hep3Vector getCrystalOriginByRO(int) const;
      // Get crystal X-axis (from front to readout side) in Mu2e coordinates
      CLHEP::Hep3Vector getCrystalAxisByRO(int) const;

    protected:

      int _nVane;
      int _nCrystalZ;
      int _nCrystalR;

      double _crystalHL;
      double _crystalHW;

      CLHEP::Hep3Vector _origin;

      double _rMin;

      double _vaneTheta;

      double _wrapperHalfThickness;
      double _roHalfTrans;
      double _roHalfThickness;

      int _nROPerCrystal;

      std::vector<Vane> _vanes;

      double _nonUniformity;
      double _timeGap;
      double _electronEdep; // energy deposition of charged particle crossing APD
      double _electronEmin; // minimum energy deposition of charged particle crossing APD

    };

} //namespace mu2e

#endif /* CalorimeterGeom_Calorimeter_hh */
