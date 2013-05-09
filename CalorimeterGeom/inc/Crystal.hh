#ifndef CalorimeterGeom_Crystal_hh
#define CalorimeterGeom_Crystal_hh
// $Id: Crystal.hh,v 1.14 2013/05/09 23:14:14 echenard Exp $
// $Author: echenard $
// $Date: 2013/05/09 23:14:14 $
//
// Hold information about position of a hexagonal cell
//
// Original author B Echenard 
//

//C++ includes
#include <vector>

//CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

     class Crystal {

	  public:

             Crystal(int id, CLHEP::Hep3Vector position) : 
	       _id(id), _position(position), _neighbours()  
	     {}


             int id(void) const                                       {return _id;}
             CLHEP::Hep3Vector const& position(void) const            {return _position;}

             std::vector<int> const& nearestNeighbours(void) const    {return _neighbours;}
             void setNearestNeighbours(std::vector<int> list)         {_neighbours = list;}


	 private:
	 
	     int                 _id;
             CLHEP::Hep3Vector   _position;
	     std::vector<int>    _neighbours;

     };

}


#endif /* CalorimeterGeom_Crystal_hh */
