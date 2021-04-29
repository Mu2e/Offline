#include "CaloCluster/inc/ClusterUtils.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "RecoDataProducts/inc/CaloHit.hh"

#include <iostream>
#include <list>


namespace mu2e {
    
    double ClusterUtils::e1() const
    {        
       return hits_[0]->energyDep();
    }

    double ClusterUtils::e2() const
    {        
       return (hits_.size()>1) ? hits_[0]->energyDep()+hits_[1]->energyDep() : hits_[0]->energyDep();
    }

    double ClusterUtils::e9() const
    {        
       const auto& neighborsId = cal_.crystal(hits_[0]->crystalID()).neighbors();

       double e9 = e1();
       for (const auto& hit : hits_)
       {
          if (std::find(neighborsId.begin(), neighborsId.end(), hit->crystalID()) != neighborsId.end()) {e9 += hit->energyDep();}
       }
       return e9;
    }

    double ClusterUtils::e25() const
    {        
       const auto& nneighborsId = cal_.crystal(hits_[0]->crystalID()).nextNeighbors();

       double e25 = e9();
       for (const auto& hit : hits_)
       {
          if (std::find(nneighborsId.begin(), nneighborsId.end(), hit->crystalID()) != nneighborsId.end()) {e25 += hit->energyDep();}
       }
       return e25;
    }




    CLHEP::Hep3Vector ClusterUtils::cog3Vector() const
    {        
        double sx(0),sy(0),sz(0),sx2(0),sy2(0),sw(0);
        fill(sx,sy,sz,sx2,sy2,sw);
        
        int iSection  = cal_.crystal(hits_[0]->crystalID()).diskID();
        CLHEP::Hep3Vector cogMu2eFrame(sx/sw,sy/sw,sz/sw);
        return cal_.geomUtil().mu2eToDiskFF(iSection,cogMu2eFrame);
    }
    
    double ClusterUtils::secondMoment() const
    {        
        double sx(0),sy(0),sz(0),sx2(0),sy2(0),sw(0);
        fill(sx,sy,sz,sx2,sy2,sw);
        return sx2-sx*sx/sw + sy2-sy*sy/sw;
    }
    
    void ClusterUtils::fill(double& sx, double& sy, double& sz, double& sx2, double& sy2, double& sw) const
    {
        int iSection  = cal_.crystal(hits_[0]->crystalID()).diskID();
        for (const auto& hit : hits_)
        {
           int    crId(hit->crystalID());
           double energy(hit->energyDep());

           if (cal_.crystal(crId).diskID() != iSection) continue;

           double xCrystal = cal_.crystal(crId).position().x();
           double yCrystal = cal_.crystal(crId).position().y();
           double zCrystal = cal_.crystal(crId).position().z();

           double weight = energy;
           if (mode_ == Logarithm)   weight = log(energy);

           sw  += weight;
           sx  += xCrystal*weight;
           sy  += yCrystal*weight;
           sz  += zCrystal*weight;
           sx2 += xCrystal*xCrystal*weight;
           sy2 += yCrystal*yCrystal*weight;
       }
    }

}


