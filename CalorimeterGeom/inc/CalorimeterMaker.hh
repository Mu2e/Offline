#ifndef CALORIMETERMAKER_HH
#define CALORIMETERMAKER_HH
// $Id: CalorimeterMaker.hh,v 1.8 2010/05/18 22:06:19 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/05/18 22:06:19 $

// original authors Julie Managan and Robert Bernstein

//
// C++ includes
#include <iomanip>
#include <memory>
#include <iostream>
#include <cmath>
#include <vector>
#include <string>

//
//Mu2e includes
#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CalorimeterGeom/inc/Crystal.hh"
#include "CalorimeterGeom/inc/CrystalDetail.hh"

#include "CalorimeterGeom/inc/Vane.hh"
#include "CalorimeterGeom/inc/RSlice.hh"


//
// other includes
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e{

  class SimpleConfig;

  namespace calorimeter{

    //
    //forward declarations
    class Calorimeter;

    class CalorimeterMaker{

    public:
      CalorimeterMaker(SimpleConfig const& config);
  
      ~CalorimeterMaker();
    

      //
      // number of vanes
      uint32_t numberOfVanes;

      //
      //half length of small side face; crystals assumed square
      double crystalHalfTrans;

      //
      //half length of long crystal axis
      double crystalHalfLong;

      //
      //number of rslices, constant distance from z-axis
      uint32_t nCrystalRSlices;

      //
      //number of zslices, constant z
      uint32_t nCrystalZSlices;

      //
      // center of calorimeter
      CLHEP::Hep3Vector calorimeterCenter;

      //
      // offset of calorimeter wrt above
      CLHEP::Hep3Vector calorimeterCenterOffset;

      // 
      // inner inscribed circle
      double rInscribed;


      // Overall azimuthal rotation
      double phi0;

      //
      // Overall tilt of calorimeter system wrt z-axis
      double theta0;



      //
      // Individual rotations in euler notation of vanes from central axis (so in system with theta_0 = 0)
      // these are supposed to be perturbations from overall rotation of vane
      std::vector<double> calorimeterVaneRotationsPhi;
      std::vector<double> calorimeterVaneRotationsTheta;
      std::vector<double> calorimeterVaneRotationsPsi;

      //
      // Individual offsets of vanes from ideal locations
      std::vector<double> calorimeterVaneOffsetsX;
      std::vector<double> calorimeterVaneOffsetsY;
      std::vector<double> calorimeterVaneOffsetsZ;

      std::string crystalMaterial;
      std::string crystalWrapper;
      double crystalWrapperHalfThickness;

      //
      // simple dumb vector
      std::vector <Crystal> allCrystals;

      // Accessor and auto_ptr to calorimeter needed by GeometryService.
      std::auto_ptr<Calorimeter> _calorimeter;
      std::auto_ptr<Calorimeter> getCalorimeterPtr() { return _calorimeter; }

    private:

      bool CheckIt();
      bool CheckVaneConsistency();

      void BuildIt();
      void MakeGenericCrystals();
      void MakeVanes();
      void MakeCalorimeter();

      //
      //some helper functions
      uint32_t NumberOfCrystalsPerVane();
      uint32_t TotalNumberOfCrystals();

      void FillNearestNeighbours();
      void FillPointersAndIndices();


      //
      // details of crystal, allowing different kinds eventually
      std::vector <CrystalDetail> _crystalDetail;


      /*
      //
      // some helper functions
      void crystalPrinter (const Crystal& s);
      void crystalPrinter2(const Crystal& s);
      void layerPrinter(const Layer& layerId);
      void devicePrinter(const Layer& deviceId);
      */
    };

  } //namespace calorimeter
} //namespace mu2e

#endif
