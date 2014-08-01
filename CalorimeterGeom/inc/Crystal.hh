#ifndef CalorimeterGeom_Crystal_hh
#define CalorimeterGeom_Crystal_hh

// $Id: Crystal.hh,v 1.15 2014/08/01 20:57:44 echenard Exp $
// $Author: echenard $
// $Date: 2014/08/01 20:57:44 $
//
// Hold information about a crystal
// Neighbors and position are given in the "mu2e" frame, not in the local disk frame
// localId and localPosition are given in the section frame (i.e. local frame)
//
// Original author B Echenard 
//



//C++ includes
#include <vector>
#include <iostream>

//CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"


namespace mu2e {

     class Crystal {

	  public:

             Crystal(int localId, CLHEP::Hep3Vector localPosition) : 
	       _localId(localId), _localPosition(localPosition),_position(), _neighborsLevel1(),_neighborsLevel2(),_neighborsLevel3()
	     {}


             int localId()                             const        {return _localId;}
             CLHEP::Hep3Vector const& localPosition()  const        {return _localPosition;}

             CLHEP::Hep3Vector const& position()       const        {return _position;}
             void setPosition(CLHEP::Hep3Vector const& pos)         {_position = pos;}


             std::vector<int> const& neighborsLevel1() const        {return _neighborsLevel1;}
             std::vector<int> const& neighborsLevel2() const        {return _neighborsLevel2;}
             std::vector<int> const& neighborsLevel3() const        {return _neighborsLevel3;}
	     
             void setNeighborsLevel1(std::vector<int> const& list)  {_neighborsLevel1 = list;}
             void setNeighborsLevel2(std::vector<int> const& list)  {_neighborsLevel2 = list;}
             void setNeighborsLevel3(std::vector<int> const& list)  {_neighborsLevel3 = list;}



	 private:
	 
	     int                 _localId;
             CLHEP::Hep3Vector   _localPosition;
             CLHEP::Hep3Vector   _position;
	     std::vector<int>    _neighborsLevel1;
	     std::vector<int>    _neighborsLevel2;
	     std::vector<int>    _neighborsLevel3;

     };

}


#endif /* CalorimeterGeom_Crystal_hh */
