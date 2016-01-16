#ifndef CalorimeterGeom_Crystal_hh
#define CalorimeterGeom_Crystal_hh

// $Id: Crystal.hh,v 1.15 2014/08/01 20:57:44 echenard Exp $
// $Author: echenard $
// $Date: 2014/08/01 20:57:44 $
//
// Hold information about a crystal
// Neighbors and position are given in the "Mu2e" frame
// localId and localPosition are given in the section frame (i.e. local frame)
//
// Original author B Echenard 
//


#include "CLHEP/Vector/ThreeVector.h"

#include <vector>
#include <iostream>


namespace mu2e {

     class Crystal {

	  public:

             Crystal(int localId, int sectionId, const CLHEP::Hep3Vector &localPosition, const CLHEP::Hep3Vector &localPositionFF) : 
	        _localId(localId), 
	        _sectionId(sectionId), 
	        _localPosition(localPosition), 
	        _localPositionFF(localPositionFF), 
	        _position(), 
	        _neighbors(), 
	        _nextNeighbors()
	     {}


             int   localId()                             const      {return _localId;}
             int   sectionId()                           const      {return _sectionId;}
             const CLHEP::Hep3Vector& localPosition()    const      {return _localPosition;}
             const CLHEP::Hep3Vector& localPositionFF()  const      {return _localPositionFF;}

             const CLHEP::Hep3Vector& position()         const      {return _position;}
             const std::vector<int>&  neighbors()        const      {return _neighbors;}	     
             const std::vector<int>&  nextNeighbors()    const      {return _nextNeighbors;}	     

             void setPosition(const CLHEP::Hep3Vector &pos)         {_position = pos;}
             void setNeighbors(const std::vector<int> &list)        {_neighbors = list;}
             void setNextNeighbors(const std::vector<int> &list)    {_nextNeighbors = list;}


	 private:
	 
	     int                 _localId;
	     int                 _sectionId;
             CLHEP::Hep3Vector   _localPosition;
             CLHEP::Hep3Vector   _localPositionFF;
             CLHEP::Hep3Vector   _position;
	     std::vector<int>    _neighbors;
	     std::vector<int>    _nextNeighbors;

     };

}


#endif 
