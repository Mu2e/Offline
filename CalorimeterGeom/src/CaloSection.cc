//
// Hold information about position of a hexagonal cell
//
// Original author B Echenard - P. Ongmongkolkul
//

// Notes: CrystalMap tesselates a plane with hexagons or squares, then selects crystals inside the ring


// C++ includes
#include <iostream>
#include <fstream>

// Mu2e includes
#include "CalorimeterGeom/inc/CaloSection.hh"

using std::cout;
using std::endl;


namespace mu2e {


      void CaloSection::setBoundsInTracker(CLHEP::Hep3Vector  const& trackerOffset, double z0, double z1, double r0, double r1) 
      {
         _zDownInTracker = z0 - trackerOffset.z();
	 _zUpInTracker   = z1 - trackerOffset.z();
	 _rInTracker     = r0;
	 _rOutTracker    = r1;
	  
      }  

      void CaloSection::print() 
      {
          cout<<"Calo section           "<<_id<<endl;
          cout<<"origin                 "<<_originLocal<<endl;
          cout<<"origin Mu2e            "<<_origin<<endl;
          cout<<"size                   "<<_size<<endl;
          cout<<"originToCrystalOrigin  "<<_originToCrystalOrigin<<endl;
          cout<<"z Down tracker         "<<_zDownInTracker<<endl;
          cout<<"z Up tracker           "<<_zUpInTracker<<endl;
          cout<<"r In tracker           "<<_rInTracker<<endl;
          cout<<"r Out tracker          "<<_rOutTracker<<endl;
      }  

 
    
}

