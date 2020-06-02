#ifndef CaloCluster_ClusterMoments_HH_
#define CaloCluster_ClusterMoments_HH_

#include "RecoDataProducts/inc/CaloCluster.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"

#include <iostream>
#include <list>

namespace mu2e {

    class ClusterMoments
    {
	public:
	    enum cogtype {Linear,LinearMod,Sqrt,Logarithm};

	    ClusterMoments(const Calorimeter& cal, const CaloCluster& caloCluster, int iSection) : 
	      cal_(cal), caloCluster_(caloCluster),iSection_(iSection),cog_(CLHEP::Hep3Vector()),secondMoment_(0.0),angle_(0.0) 
	    {};

	    void calculate(cogtype mode = Linear); 

            const CLHEP::Hep3Vector&   cog         () const {return cog_;}
                  double               secondMoment() const {return secondMoment_;}
                  double               angle       () const {return angle_;}



	private:
            const Calorimeter&  cal_;
	    const CaloCluster&  caloCluster_;
	    int                 iSection_;	     
	    CLHEP::Hep3Vector   cog_;
	    double              secondMoment_;
	    double              angle_;
    };


} 

#endif
