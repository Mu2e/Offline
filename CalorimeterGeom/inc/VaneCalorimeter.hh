#ifndef CalorimeterGeom_VaneCalorimeter_hh
#define CalorimeterGeom_VaneCalorimeter_hh
//
// Hold all geometry and identifier information about
// a Calorimeter.  In order to insulate this class from
// knowledge of databases etc, this class must not know
// how to make itself.
//
// $Id: VaneCalorimeter.hh,v 1.1 2012/09/08 02:24:25 echenard Exp $
// $Author: echenard $
// $Date: 2012/09/08 02:24:25 $
//
// Original author R. Bernstein and Rob Kutschke
//

#include <vector>
#include <iostream>

//
// Mu2e includes
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CalorimeterGeom/inc/Vane.hh"

namespace mu2e {
    
    class VaneCalorimeter: public Calorimeter{

      friend class VaneCalorimeterMaker;

    public:
      VaneCalorimeter(){}
      ~VaneCalorimeter(){}

      CLHEP::Hep3Vector const& getOrigin(void) const { return _origin; }
      double innerRaidus(void) const                 { return _rMin; }
      double outherRadius(void) const                { return _rMax;}


      unsigned int nVane(void) const           { return _nVane; }
      Vane const& getVane(int i) const         { return _vanes.at(i); }

      unsigned int nCrystalPerVane(void) const { return _nCrystalZ*_nCrystalR; }
      unsigned int nCrystal(void) const        { return _nVane*_nCrystalZ*_nCrystalR; }
      unsigned int nCrystalR(void) const       { return _nCrystalR; }
      unsigned int nCrystalZ(void) const       { return _nCrystalZ; }

      unsigned int nRO(void) const             { return _nROPerCrystal*_nVane*_nCrystalZ*_nCrystalR; }
      unsigned int nROPerCrystal(void) const   { return _nROPerCrystal; }

      double crystalHalfSize(void) const       { return _crystalHW; }
      double crystalHalfLength(void) const     { return _crystalHL; }
      double wrapperThickness(void) const      { return _wrapperThickness; }
      double shellThickness(void) const        { return _shellThickness;}
      double roHalfSize(void) const            { return _roHalfTrans; }
      double roHalfThickness(void) const       { return _roHalfThickness; }

      double getNonuniformity(void) const      { return _nonUniformity; }
      double getTimeGap(void) const            { return _timeGap; }
      double getElectronEdep(void) const       { return _electronEdep; }
      double getElectronEmin(void) const       { return _electronEmin; }
      
      double getApdMeanNoise(void) const       { return _apdMeanNoise;}
      double getApdSigmaNoise(void) const      { return _apdSigmaNoise;}
      
      double getLysoLightYield(void) const     {return _lysoLightYield;}
      double getApdQuantumEff(void) const      {return _apdQuantumEff;}
      double getAPDCollectEff(void) const      {return _lightCollectEffAPD;}

      // Vane ID (0..nVanes-1)
      int getVaneByRO(int roid) const {
        return roid/(_nCrystalZ*_nCrystalR*_nROPerCrystal);
      }
      
      // Crystal ID (0..Number_of_crystals_in_calorimeter-1)
      int getCrystalByRO(int roid) const { return (roid/_nROPerCrystal); }
      
      //RO base id = crystal_id*nROPerCryastal ... crystal_id*nROPerCryastal+ nROPerCryastal-1
      int getROBaseByCrystal(int id) const {return (id*_nROPerCrystal);}
      
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

      // Transfer Mu2e coordinates to local vane coordinates
      CLHEP::Hep3Vector toVaneFrame(int, CLHEP::Hep3Vector const&) const;

      //inverse transformation of "toVaneFrame"
      CLHEP::Hep3Vector fromVaneFrame(int vaneid, CLHEP::Hep3Vector const& pos) const;

      // Get crystal center in Mu2e coordinates
      CLHEP::Hep3Vector getCrystalOriginByRO(int) const;

      // Get crystal X-axis (from front to readout side) in Mu2e coordinates
      CLHEP::Hep3Vector getCrystalAxisByRO(int) const;



    protected:

      std::vector<Vane> _vanes;
      int    _nVane;
      int    _nCrystalZ;
      int    _nCrystalR;
      double _rMin;
      double _rMax;
      double _vaneTheta;

      CLHEP::Hep3Vector _origin;

      double _crystalHL;
      double _crystalHW;

      double _wrapperThickness;
      double _shellThickness;

      int    _nROPerCrystal;
      double _roHalfTrans;
      double _roHalfThickness;

      double _nonUniformity;
      double _timeGap;
      double _electronEdep; // energy deposition of charged particle crossing APD
      double _electronEmin; // minimum energy deposition of charged particle crossing APD

      double _apdMeanNoise; //MeV
      double _apdSigmaNoise;//MeV

      double _lysoLightYield;
      double _apdQuantumEff;//quantum efficiency for Hamamatsu S8664-1010 for a radiation wavelenght of 402nm (typical of lyso)
      double _lightCollectEffAPD;//light collection efficiency for 30 × 30 mm2 area of the crystal efficiency for Hamamatsu S8664-1010

    };

} //namespace mu2e

#endif /* CalorimeterGeom_VaneCalorimeter_hh */
