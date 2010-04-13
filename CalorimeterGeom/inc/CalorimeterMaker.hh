#ifndef CALORIMETERMAKER_HH
#define CALORIMETERMAKER_HH
// $Id: CalorimeterMaker.hh,v 1.3 2010/04/13 17:14:41 rhbob Exp $
// $Author: rhbob $
// $Date: 2010/04/13 17:14:41 $

// original authors Julie Managan and Robert Bernstein
namespace mu2e{
  namespace calorimeter{

    //
    // C++ includes
#include <vector>

    //
    //Mu2e includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CalorimeterGeom/inc/Device.hh"
#include "CalorimeterGeom/inc/LayerInfo.hh"
#include "CalorimeterGeom/inc/Layer.hh"


    //
    //forward declarations
    class Calorimeter;

    class CalorimeterMaker{

    public:
      CalorimeterMaker( int nVanes,
			std::vector<LayerInfo> vaneInfo,
			double r0,
			double rInscribed,
			double caloZLength,
			double halfSide,
			CLHEP::Hep3Vector center,
			double phi0
			);
  
      ~CalorimeterMaker ();
  
      const Calorimeter& getCalorimeter() const { return *_calo;}
  
    private:

      void CheckFit();
      void CheckVaneConsistency();
      int  totalCrystals() const;

      void MakeVanes();
      void MakeDetails();
      void FillNearestNeighbours();

      // Number of vanes in the calorimeter
      int _nVanes;

      // Number of crystals per layer in the vanes.
      std::vector<LayerInfo> _vaneInfo;

      // Radius to the middle of each vane.
      double _r0;
      double _rInscribed;

      // Length of the calorimeter in Z
      double _caloZLength;
  
      // Crystal halfSide.
      double _halfSide;

      // Center of the calorimeter
      CLHEP::Hep3Vector _center;

      // Overall azimuthal rotation
      double _phi0;

      // The actual Calorimeter.
      Calorimeter* _calo;

    };
  } //namespace calorimeter
} //namespace mu2e

#endif
