#include "CaloMC/inc/ShowerStepUtil.hh"
#include "cetlib_except/exception.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/SymMatrix.h"

#include <vector>
#include <iostream>
#include <algorithm>


namespace mu2e {
    
    void ShowerStepUtil::add(unsigned i, double eDepG4, double eDepVis, double time, double momentum, CLHEP::Hep3Vector& pos)
    {               
        if (i > imax_) throw cet::exception("Rethrow")<< "[CaloMC/ShowerStepUtil] Index out of bound " << i << std::endl;

        //init buffer if needed
        if (n_[i]==0) {pIn_[i] = momentum; t0_[i] = time;} 
        
        double weight = (type_ == weight_type::energy)  ? eDepG4 : 1.0;
        
        n_[i]      += 1;
        eDepG4_[i] += eDepG4;
        eDepVis_[i]+= eDepVis;
        pIn_[i]     = std::max(pIn_[i],momentum);
        time_[i]   += time*weight;             
        x_[i]      += pos.x()*weight;             
        y_[i]      += pos.y()*weight;             
        z_[i]      += pos.z()*weight;             
        w_[i]      += weight;             
    }


    void ShowerStepUtil::reset(unsigned i)
    {
        if (i > imax_) throw cet::exception("Rethrow")<< "[CaloMC/ShowerStepUtil] Index out of bound " << i << std::endl;      
        n_[i]       = 0;
        eDepG4_[i]  = 0;
        eDepVis_[i] = 0;
        pIn_[i]     = 0;
        time_[i]    = x_[i] = y_[i] = z_[i] = w_[i] = 0;
    }

    
    CLHEP::Hep3Vector& ShowerStepUtil::pos(unsigned i)
    {        
        if (i > imax_) throw cet::exception("Rethrow")<< "[CaloMC/ShowerStepUtil] Index out of bound " << i << std::endl;
        
        pos_[0] = x_[i]/w_[i];
        pos_[1] = y_[i]/w_[i];
        pos_[2] = z_[i]/w_[i];
        return pos_;
    }


    void ShowerStepUtil::printBucket(unsigned i)
    {
        if (i > imax_) throw cet::exception("Rethrow")<< "[CaloMC/ShowerStepUtil] Index out of bound " << i << std::endl;        
        std::cout<<"Entries= "<<n_[i]<<" Energy = "<<eDepG4_[i]<<" Time = "<<time_[i]/w_[i]
                 <<" pos=("<<x_[i]/w_[i]<<","<<y_[i]/w_[i]<<","<<z_[i]/w_[i]<<")  momentum="<<pIn_[i]<<std::endl;
    }

}



