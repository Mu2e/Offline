#ifndef CalorimeterGeom_CalorimeterMaker_hh
#define CalorimeterGeom_VaneCalorimeterMaker_hh
// $Id: VaneCalorimeterMaker.hh,v 1.4 2013/03/15 15:52:03 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/03/15 15:52:03 $

// original authors Julie Managan and Robert Bernstein

// C++ includes
#include <iomanip>
#include <memory>
#include <iostream>
#include <cmath>
#include <vector>
#include <string>

//Mu2e includes
#include "CalorimeterGeom/inc/VaneCalorimeter.hh"
#include "CalorimeterGeom/inc/Vane.hh"
#include "ConfigTools/inc/SimpleConfig.hh"

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"



namespace mu2e{


  class SimpleConfig;

    
    //forward declarations
    class VaneCalorimeter;

    class VaneCalorimeterMaker{

    public:
     
      VaneCalorimeterMaker(SimpleConfig const& config, double solenoidOffset);
      ~VaneCalorimeterMaker();
      
      
      // Accessor and unique_ptr to VaneCalorimeter needed by GeometryService.
      std::unique_ptr<VaneCalorimeter> _calo;
      std::unique_ptr<VaneCalorimeter> calorimeterPtr() { return std::move(_calo); }

    private:

       void CheckIt(void);
       void MakeVanes(void);

       int _verbosityLevel;
   };

} 

#endif 
