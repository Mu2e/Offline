//
// Original author B.Echenard
//

#ifndef CaloCluster_ClusterMoments_HH_
#define CaloCluster_ClusterMoments_HH_


// Mu2e includes
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"


// C++ includes
#include <iostream>
#include <list>

using CLHEP::Hep3Vector;


namespace mu2e {


    class ClusterMoments {


	 public:
             
	     enum cogtype {Linear,LinearMod,Sqrt,Logarithm};
             
	     ClusterMoments(Calorimeter const& cal, CaloCluster const& caloCluster, int iSection) : 
	       cal_(cal), caloCluster_(caloCluster),iSection_(iSection),cog_(CLHEP::Hep3Vector(0,0,0)),secondMoment_(0),angle_(0) 
	     {};

             ~ClusterMoments(){};

	                  
	     void calculate(cogtype mode = Linear); 


             CLHEP::Hep3Vector const& cog()     const {return cog_;}
             double            secondMoment()   const {return secondMoment_;}
             double            angle()          const {return angle_;}



	 private:

             Calorimeter const&  cal_;
	     CaloCluster const&  caloCluster_;
	     int                 iSection_;
	     
	     CLHEP::Hep3Vector   cog_;
	     double              secondMoment_;
	     double              angle_;
	     
	     
	     
    };


} 

#endif
