#include "CaloMC/inc/ShowerStepUtil.hh"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include <vector>
#include <iostream>
#include <algorithm>


namespace mu2e {
    
    
    void ShowerStepUtil::init(int i, double time, double momentum, const CLHEP::Hep3Vector& posIn)
    { 
        t0_.at(i) = time; 
        pIn_.at(i) = momentum; 
        posIn_.at(i) = posIn;
    }


    void ShowerStepUtil::add(int i, double edep, double time, double momentum, CLHEP::Hep3Vector& pos)
    {               
        double weight = (type_ == weight_type::energy)  ? edep : 1.0;

        n_.at(i)    += 1;
        edep_.at(i) += edep;
        pIn_.at(i)  = std::max(pIn_.at(i),momentum);

        time_.at(i) += time*weight;             
        x_.at(i)    += pos.x()*weight;             
        y_.at(i)    += pos.y()*weight;             
        z_.at(i)    += pos.z()*weight;             

        x2_.at(i)   += pos.x()*pos.x()*weight;             
        y2_.at(i)   += pos.y()*pos.y()*weight;             
        z2_.at(i)   += pos.z()*pos.z()*weight;             

        xy_.at(i)   += pos.x()*pos.y()*weight;             
        xz_.at(i)   += pos.x()*pos.z()*weight;             
        yz_.at(i)   += pos.y()*pos.z()*weight;             

        w_.at(i)    += weight;             
        w2_.at(i)   += weight*weight;             
    }


    void ShowerStepUtil::reset(int i)
    {
        n_[i]     = 0;
        edep_[i]  = 0;
        time_[i]  = 0;
        pIn_[i]   = 0;
        x_[i]  = y_[i]  = z_[i]  = 0;
        x2_[i] = y2_[i] = z2_[i] = 0;
        xy_[i] = yz_[i] = xz_[i] = 0;
        w_[i]  = w2_[i] = 0;        
    }

    
    CLHEP::Hep3Vector& ShowerStepUtil::pos(int i)
    {        
        pos_[0] = x_[i]/w_[i];
        pos_[1] = y_[i]/w_[i];
        pos_[2] = z_[i]/w_[i];
      
        return pos_;
    }

    CLHEP::HepSymMatrix& ShowerStepUtil::covPos(int i)
    {        
        double norm = w_[i]/(w_[i]*w_[i]-w2_[i]);
      
        cov_[0][0] = norm*(x2_[i] - x_[i]*x_[i]/w_[i]);
        cov_[1][1] = norm*(y2_[i] - y_[i]*y_[i]/w_[i]);
        cov_[2][2] = norm*(z2_[i] - z_[i]*z_[i]/w_[i]);
        cov_[0][1] = norm*(xy_[i] - x_[i]*y_[i]/w_[i]);
        cov_[0][2] = norm*(xz_[i] - x_[i]*z_[i]/w_[i]);
        cov_[1][2] = norm*(yz_[i] - y_[i]*z_[i]/w_[i]);

        return cov_;
    }


    void ShowerStepUtil::printBucket(int i)
    {
        std::cout<<"Entries= "<<n_[i]<<" Energy = "<<edep_[i]<<" Time = "<<time_[i]/w_[i]<<" pos=("<<x_[i]/w_[i]<<","<<y_[i]/w_[i]<<","<<z_[i]/w_[i]
                 <<")  momentum="<<pIn_[i]<<std::endl;
    }

}



