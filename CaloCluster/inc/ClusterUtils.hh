#ifndef CaloCluster_ClusterUtils_HH_
#define CaloCluster_ClusterUtils_HH_

#include "RecoDataProducts/inc/CaloCluster.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"

#include <iostream>
#include <list>

namespace mu2e {

    class ClusterUtils
    {
        public:
            enum Cogtype {Linear, Logarithm};

            ClusterUtils(const Calorimeter& cal, const CaloCluster& cluster, Cogtype mode = Linear) : 
                cal_(cal), hits_(cluster.caloHitsPtrVector()), mode_(mode)
            {};

            ClusterUtils(const Calorimeter& cal, const CaloHitPtrVector& hits, Cogtype mode = Linear) : 
                cal_(cal), hits_(hits), mode_(mode)
            {};

            CLHEP::Hep3Vector cog3Vector  () const;
            double            secondMoment() const;
            double            e1          () const;
            double            e2          () const;
            double            e9          () const;
            double            e25         () const;


        private:       
           const Calorimeter&     cal_;
           const CaloHitPtrVector hits_;
           const Cogtype          mode_;
           
           void fill(double& sx, double& sy, double& sz, double& sx2, double& sy2, double& sw) const;
   };
} 

#endif
