//
// Utility to compute the center-of-gravity and other dynamical cluster quantities
//  
// Original author B. Echenard
//

// Mu2e includes
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "CaloCluster/inc/CaloClusterCogCalculator.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"


// C++ includes
#include <iostream>
#include <list>

using CLHEP::Hep3Vector;


namespace mu2e {


       CaloClusterCogCalculator::CaloClusterCogCalculator() : _offset(4) {}


       CLHEP::Hep3Vector CaloClusterCogCalculator::calculateCog(std::list<CaloCrystalHit const*> cluster, cogtype mode) 
       {
           
	   Calorimeter const & cal = *(GeomHandle<Calorimeter>());
	   
	   CLHEP::Hep3Vector aVector(0,0,0);
	   double sumWeights(0);    

	   for (std::list<CaloCrystalHit const*>::const_iterator it = cluster.begin(); it != cluster.end(); ++it)
	   {        
	     double energy = (*it)->energyDep();
	     CLHEP::Hep3Vector crystalPos = cal.crystalOrigin((*it)->id());

	     if (energy > 1e-6) {
        	 double weight = energy;
		 if (mode == logartihm) weight = _offset+log(energy);
		 aVector[0] += crystalPos.x()*weight;
		 aVector[1] += crystalPos.y()*weight;
		 aVector[2] += crystalPos.z()*weight;
		 sumWeights += weight;
	     }

	  }

	  aVector[0] /= sumWeights;
	  aVector[1] /= sumWeights;
	  aVector[2] /= sumWeights;

	  return aVector;

       }


}
