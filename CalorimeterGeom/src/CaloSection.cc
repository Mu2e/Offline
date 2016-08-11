//
// Hold information about position of a hexagonal cell
//
// Original author B Echenard - P. Ongmongkolkul
//

// Notes: CrystalMap tesselates a plane with hexagons or squares, then selects crystals inside the ring


// C++ includes
#include <iostream>

// Mu2e includes
#include "CalorimeterGeom/inc/CaloSection.hh"



namespace mu2e {



      void CaloSection::print(std::ostream &os) const
      {
           os<<"Calo section           "<<_id<<std::endl;
           os<<"origin                 "<<_originLocal<<std::endl;
           os<<"origin Mu2e            "<<_origin<<std::endl;
           os<<"size                   "<<_size<<std::endl;
           os<<"rotation               "<<_rotation<<std::endl;
           os<<"originToCrystalOrigin  "<<_originToCrystalOrigin<<std::endl;
           os<<"z Front                "<<_frontFaceCenter.z()<<std::endl;
           os<<"z Back                 "<<_backFaceCenter.z()<<std::endl;
           os<<"r In tracker           "<<_innerRadius<<std::endl;
           os<<"r Out tracker          "<<_outerRadius<<std::endl;
      }  

 
    
}

