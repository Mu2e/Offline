#ifndef CalorimeterGeom_VaneCalorimeter_hh
#define CalorimeterGeom_VaneCalorimeter_hh
//
// Hold all geometry and identifier information about
// a Calorimeter.  In order to insulate this class from
// knowledge of databases etc, this class must not know
// how to make itself.
//
// $Id: VaneCalorimeter.hh,v 1.2 2012/11/17 00:06:25 echenard Exp $
// $Author: echenard $
// $Date: 2012/11/17 00:06:25 $
//
// Original author R. Bernstein and Rob Kutschke
//

//C++ includes
#include <vector>

// Mu2e includes
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CalorimeterGeom/inc/Vane.hh"

//CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"
    

namespace mu2e {

    class VaneCalorimeter: public Calorimeter{

      friend class VaneCalorimeterMaker;


    public:

      VaneCalorimeter(){}
      ~VaneCalorimeter(){}


      int nVane(void) const                    {return _nVane; }
      Vane const& getVane(int i) const         {return _vanes.at(i); }

      int nCrystalR(void) const                {return _nCrystalR; }
      int nCrystalZ(void) const                {return _nCrystalZ; }

      unsigned int nRO(void) const             {return _nROPerCrystal*_nVane*_nCrystalZ*_nCrystalR; }
      unsigned int nROPerCrystal(void) const   {return _nROPerCrystal; }
      unsigned int nCrystal(void) const        {return _nVane*_nCrystalZ*_nCrystalR; }

      double crystalHalfSize(void) const       {return _crystalHW; }
      double crystalHalfLength(void) const     {return _crystalHL; }
      double crystalVolume(void) const         {return 8*_crystalHW*_crystalHL;}

      double wrapperThickness(void) const      {return _wrapperThickness; }
      double shellThickness(void) const        {return _shellThickness;}
      double roHalfSize(void) const            {return _roHalfTrans; }
      double roHalfThickness(void) const       {return _roHalfThickness; }

      double getNonuniformity(void) const      {return _nonUniformity; }
      double getTimeGap(void) const            {return _timeGap; }
      double getElectronEdep(void) const       {return _electronEdep; }
      double getElectronEmin(void) const       {return _electronEmin; }
      
      double getApdMeanNoise(void) const       {return _apdMeanNoise;}
      double getApdSigmaNoise(void) const      {return _apdSigmaNoise;}
      
      double getLysoLightYield(void) const     {return _lysoLightYield;}
      double getApdQuantumEff(void) const      {return _apdQuantumEff;}
      double getAPDCollectEff(void) const      {return _lightCollectEffAPD;}




      CLHEP::Hep3Vector const& getOrigin(void) const { return _origin; }      
      CLHEP::Hep3Vector toCrystalFrame(int crystalId,   CLHEP::Hep3Vector const& pos) const; //transform Mu2e -> local coordinates
      CLHEP::Hep3Vector toSectionFrame(int sectionId,   CLHEP::Hep3Vector const& pos) const;  //transform Mu2e -> local  coordinates
      CLHEP::Hep3Vector fromSectionFrame(int sectionId, CLHEP::Hep3Vector const& pos) const; //transform local -> Mu2e coordinates
      CLHEP::Hep3Vector getCrystalOrigin(int crystalId) const;   //crystal center in Mu2e coordinates      
      CLHEP::Hep3Vector getLocalCrystalOrigin(int crystalId) const; //crystal center in vane frame
      CLHEP::Hep3Vector getCrystalAxis(int crystalId) const;   //crystal axis (front -> readout) in Mu2e coordinates



      // Navigating crystal and vane Id's
      int getCrystalByRO(int roid) const              {return (roid/_nROPerCrystal); }
      int getROBaseByCrystal(int crystalId) const     {return (crystalId*_nROPerCrystal);}

      int getCaloSectionId(int crystalId) const       {return crystalId/(_nCrystalZ*_nCrystalR);}
      int getLocalCrystalId(int crystalId) const      {return crystalId%(_nCrystalZ*_nCrystalR);}

      std::vector<int> getNeighbors(int crystalId, int level=1) const;




//keep only for backward compatibility, will disappear in the future.
//should disappear soon

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

      // Get crystal center in Mu2e coordinates
      CLHEP::Hep3Vector getCrystalOriginByRO(int roid) const { return getCrystalOrigin(getCrystalByRO(roid));  }

      int getVaneByRO(int roid) const {
        return roid/(_nCrystalZ*_nCrystalR*_nROPerCrystal);
      }

      // Get crystal center in Mu2e coordinates
      CLHEP::Hep3Vector getCrystalAxisByRO(int roid) const { return getCrystalAxis(getCrystalByRO(roid));  }


CLHEP::Hep3Vector toVaneFrame(int vaneId, CLHEP::Hep3Vector const& pos) const   {return toSectionFrame(vaneId,pos);}
CLHEP::Hep3Vector fromVaneFrame(int vaneId, CLHEP::Hep3Vector const& pos) const {return fromSectionFrame(vaneId,pos);}


    protected:

      std::vector<Vane> _vanes;
      int    _nVane;
      int    _nCrystalZ;
      int    _nCrystalR;
      double _rMin;

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

}

#endif /* CalorimeterGeom_VaneCalorimeter_hh */
