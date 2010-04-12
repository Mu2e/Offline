#ifndef CALORIMETERMAKER_HH
#define CALORIMETERMAKER_HH

//
// C++ includes
#include <vector>

//
//Mu2e includes
#include "CLHEP/Vector/ThreeVector.h"
#include "Calorimeter/inc/Device.hh"
#include "Calorimeter/inc/LayerInfo.hh"
#include "Calorimeter/inc/Layer.hh"


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

#endif
