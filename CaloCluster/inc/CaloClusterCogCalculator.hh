//
// Original author B.Echenard
//

#ifndef CaloCluster_CaloClusterCogCalculator_HH_
#define CaloCluster_CaloClusterCogCalculator_HH_


// Mu2e includes
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"


// C++ includes
#include <iostream>
#include <list>

using CLHEP::Hep3Vector;


namespace mu2e {


    class CaloClusterCogCalculator{

	 public:
             
	     enum cogtype {linear,sqrt,logartihm};
             
	     CaloClusterCogCalculator();
             ~CaloClusterCogCalculator(){};
	                  
	     CLHEP::Hep3Vector calculateCog(std::list<CaloCrystalHit const*> cluster, cogtype mode = logartihm); 


	 private:

             double _offset;
    };


} // end namespace mu2e

#endif
