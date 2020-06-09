#include "CaloCluster/inc/ClusterMoments.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "RecoDataProducts/inc/CaloCrystalHit.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"

#include <iostream>
#include <list>

namespace mu2e {


    void ClusterMoments::calculate(cogtype mode)
    {
        //calculate first and second moments in one pass
        double sxi(0),sxi2(0),syi(0),syi2(0),sxyi(0),szi(0),szi2(0),swi(0);
        for (const auto& hit : caloCluster_.caloCrystalHitsPtrVector())
        {
            int    crId(hit->id());
            double energy(hit->energyDep());

            if (cal_.crystal(crId).diskId() != iSection_) continue;

            double xCrystal = cal_.crystal(crId).position().x();
            double yCrystal = cal_.crystal(crId).position().y();
            double zCrystal = cal_.crystal(crId).position().z();

            double weight = energy;
            if (mode == Logarithm)   weight = log(energy);

            sxi  += xCrystal*weight;
            sxi2 += xCrystal*xCrystal*weight;
            syi  += yCrystal*weight;
            syi2 += yCrystal*yCrystal*weight;
            sxyi += xCrystal*yCrystal*weight;
            szi  += zCrystal*weight;
            szi2 += zCrystal*yCrystal*weight;
            swi  += weight;
       }


       CLHEP::Hep3Vector cogMu2eFrame(sxi/swi,syi/swi,szi/swi);
       cog_ = cal_.geomUtil().mu2eToDiskFF(iSection_,cogMu2eFrame);

       secondMoment_ = sxi2 -sxi*sxi/swi + syi2 -syi*syi/swi;

       double beta  = (sxyi - sxi*syi/swi)/(sxi2-sxi*sxi/swi);
       angle_ = atan(beta);
       if (angle_ < 0) angle_ += CLHEP::pi;

    }

}


