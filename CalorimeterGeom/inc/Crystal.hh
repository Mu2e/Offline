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

             Crystal(int localId, int sectionId, CLHEP::Hep3Vector localPosition) : 
	       _localId(localId), _sectionId(sectionId), _localPosition(localPosition), _position(), _neighbors(), _nextNeighbors()
	     {}


             int localId()                             const        {return _localId;}
             int sectionId()                           const        {return _sectionId;}
             CLHEP::Hep3Vector const& localPosition()  const        {return _localPosition;}

             CLHEP::Hep3Vector const& position()       const        {return _position;}
             std::vector<int> const& neighbors()       const        {return _neighbors;}	     
             std::vector<int> const& nextNeighbors()   const        {return _nextNeighbors;}	     

             void setPosition(CLHEP::Hep3Vector const& pos)         {_position = pos;}
             void setNeighbors(std::vector<int> const& list)        {_neighbors = list;}
             void setNextNeighbors(std::vector<int> const& list)    {_nextNeighbors = list;}



	 private:
	 
	     int                 _localId;
	     int                 _sectionId;
             CLHEP::Hep3Vector   _localPosition;
             CLHEP::Hep3Vector   _position;
	     std::vector<int>    _neighbors;
	     std::vector<int>    _nextNeighbors;

     };

}


#endif 
