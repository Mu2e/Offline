#ifndef CalorimeterGeom_VaneCalorimeter_hh
#define CalorimeterGeom_VaneCalorimeter_hh
//
// Hold all geometry and identifier information about
// a Calorimeter.  In order to insulate this class from
// knowledge of databases etc, this class must not know
// how to make itself.
//
// $Id: VaneCalorimeter.hh,v 1.3 2013/03/05 20:33:25 aluca Exp $
// $Author: aluca $
// $Date: 2013/03/05 20:33:25 $
//
// Original author R. Bernstein and Rob Kutschke
//

//C++ includes
#include <vector>

// Mu2e includes
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CalorimeterGeom/inc/Vane.hh"

// Framework includes
#include "cetlib/pow.h"

//CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"


namespace mu2e {

class VaneCalorimeter: public Calorimeter{

        friend class VaneCalorimeterMaker;


public:

        VaneCalorimeter(){}
        ~VaneCalorimeter(){}

        double innerRaidus(void) const                 { return _rMin; }
        double outherRadius(void) const                { return _rMax;}


        int nVane(void) const                    {return _nVane; }
        Vane const& vane(int i) const         {return _vanes.at(i); }

        int nCrystalR(void) const                {return _nCrystalR; }
        int nCrystalZ(void) const                {return _nCrystalZ; }

        unsigned int nRO(void) const             {return _nROPerCrystal*_nVane*_nCrystalZ*_nCrystalR; }
        unsigned int nROPerCrystal(void) const   {return _nROPerCrystal; }
        unsigned int nCrystal(void) const        {return _nVane*_nCrystalZ*_nCrystalR; }

        double crystalHalfSize(void) const       {return _crystalHW; }
        double crystalHalfLength(void) const     {return _crystalHL; }
        double crystalVolume(void) const         {return 8*_crystalHW*_crystalHL;}
        
        double shieldHalfThickness(void) const           {return _shieldHalfThickness;}
        double neutronAbsorberHalfThickness(void) const  {return _neutronAbsorberHalfThickness;}

        double envelopeRmin(void) const                  {return _vanes.at(0).innerRadius() - 1.0 ;}
        double envelopeRmax(void) const                  {return sqrt(cet::square(_vanes.at(0).outerRadius() + 1.0) + cet::square(_crystalHL) ) + 1.0;}
        double envelopeHalfLength(void) const            {return _vanes.at(0).size().z() +2.0*(_shieldHalfThickness + _neutronAbsorberHalfThickness) + 1.0;}

        
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




        CLHEP::Hep3Vector const& origin(void) const { return _origin; }
        CLHEP::Hep3Vector toCrystalFrame(int crystalId,   CLHEP::Hep3Vector const& pos) const; //transform Mu2e -> local coordinates
        CLHEP::Hep3Vector toSectionFrame(int sectionId,   CLHEP::Hep3Vector const& pos) const;  //transform Mu2e -> local  coordinates
        CLHEP::Hep3Vector fromSectionFrame(int sectionId, CLHEP::Hep3Vector const& pos) const; //transform local -> Mu2e coordinates
        CLHEP::Hep3Vector crystalOrigin(int crystalId) const;   //crystal center in Mu2e coordinates
        CLHEP::Hep3Vector localCrystalOrigin(int crystalId) const; //crystal center in vane frame
        CLHEP::Hep3Vector crystalAxis(int crystalId) const;   //crystal axis (front -> readout) in Mu2e coordinates



        // Navigating crystal and vane Id's
        int crystalByRO(int roid) const              {return (roid/_nROPerCrystal); }
        int ROBaseByCrystal(int crystalId) const     {return (crystalId*_nROPerCrystal);}

        int caloSectionId(int crystalId) const       {return crystalId/(_nCrystalZ*_nCrystalR);}
        int localCrystalId(int crystalId) const      {return crystalId%(_nCrystalZ*_nCrystalR);}

        std::vector<int> neighbors(int crystalId, int level=1) const;




        //keep only for backward compatibility, will disappear in the future.
        //should disappear soon

        // Crystal ID within a vane (0..Number_of_crystals_in_vane-1)
        int crystalVaneByRO(int roid) const {
                return (roid/_nROPerCrystal)%(_nCrystalZ*_nCrystalR);
        }

        // Crystal R-coordinate within a vane (0..nCrystalR-1)
        int crystalRByRO(int roid) const {
                return ((roid/_nROPerCrystal)%(_nCrystalZ*_nCrystalR))/_nCrystalZ;
        }

        // Crystal Z-coordinate within a vane (0..nCrystalZ-1)
        int crystalZByRO(int roid) const {
                return ((roid/_nROPerCrystal)%(_nCrystalZ*_nCrystalR))%_nCrystalZ;
        }

        // Get crystal center in Mu2e coordinates
        CLHEP::Hep3Vector crystalOriginByRO(int roid) const { return crystalOrigin(crystalByRO(roid));  }

        int vaneByRO(int roid) const {
                return roid/(_nCrystalZ*_nCrystalR*_nROPerCrystal);
        }

        // Get crystal center in Mu2e coordinates
        CLHEP::Hep3Vector crystalAxisByRO(int roid) const { return crystalAxis(crystalByRO(roid));  }


        CLHEP::Hep3Vector toVaneFrame(int vaneId, CLHEP::Hep3Vector const& pos) const   {return toSectionFrame(vaneId,pos);}
        CLHEP::Hep3Vector fromVaneFrame(int vaneId, CLHEP::Hep3Vector const& pos) const {return fromSectionFrame(vaneId,pos);}


protected:

        std::vector<Vane> _vanes;
        int    _nVane;
        int    _nCrystalZ;
        int    _nCrystalR;
        double _rMin;
        double _rMax;

        CLHEP::Hep3Vector _origin;

        double _crystalHL;
        double _crystalHW;
        
        double _shieldHalfThickness;
        double _neutronAbsorberHalfThickness;
        
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
