#ifndef CALORIMETERMAKER_HH
#define CALORIMETERMAKER_HH
// $Id: CalorimeterMaker.hh,v 1.5 2010/04/27 18:46:35 rhbob Exp $
// $Author: rhbob $
// $Date: 2010/04/27 18:46:35 $

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
/*
#include "CalorimeterGeom/inc/Device.hh"
#include "CalorimeterGeom/inc/LayerInfo.hh"
#include "CalorimeterGeom/inc/Layer.hh"
*/

//
// other includes
#include "CLHEP/Vector/ThreeVector.h"

using namespace std;

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
      //number of rows, rows defined as starting from center and going out
      uint32_t nCrystalRows;

      //
      //number of columns, orthogonal to rows
      uint32_t nCrystalColumns;

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
      // Individual rotations in x,y,z of vanes from central axis (so in system with theta_0 = 0)
      std::vector<double> calorimeterVaneRotationsX;
      std::vector<double> calorimeterVaneRotationsY;
      std::vector<double> calorimeterVaneRotationsZ;

      //
      // Individual offsets of vanes from ideal locations
      std::vector<double> calorimeterVaneOffsetsX;
      std::vector<double> calorimeterVaneOffsetsY;
      std::vector<double> calorimeterVaneOffsetsZ;

      std::string crystalMaterial;
      std::string crystalWrapper;
      double crystalWrapperThickness;

      //
      // the calorimeter and a vector of crystals in it, for dumb access
      Calorimeter* _calorimeter;
      std::vector <Crystal> allCrystals;


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
