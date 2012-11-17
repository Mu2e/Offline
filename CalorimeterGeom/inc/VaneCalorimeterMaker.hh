#ifndef CalorimeterGeom_CalorimeterMaker_hh
#define CalorimeterGeom_VaneCalorimeterMaker_hh
// $Id: VaneCalorimeterMaker.hh,v 1.2 2012/11/17 00:06:25 echenard Exp $
// $Author: echenard $
// $Date: 2012/11/17 00:06:25 $

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
      
      
      // Accessor and auto_ptr to VaneCalorimeter needed by GeometryService.
      std::auto_ptr<VaneCalorimeter> _calo;
      std::auto_ptr<VaneCalorimeter> getCalorimeterPtr() { return _calo; }

    private:

       void CheckIt(void);
       void MakeVanes(void);

    };

} //namespace mu2e

#endif /* CalorimeterGeom_VaneCalorimeterMaker_hh */
