#include "HitMakers/inc/CaloShowerStepMCUtil.hh"

#include <vector>
#include <iostream>
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/SymMatrix.h"


namespace mu2e {


    void CaloShowerStepMCUtil::add(int i, double edep, double time, CLHEP::Hep3Vector& pos)
    {       	
	 double weight = (_type == weight_type::energy)  ? edep : 1.0;

	 _n.at(i)    += 1;
	 _edep.at(i) += edep;

	 _time.at(i) += time*weight;             
	 _x.at(i)    += pos.x()*weight;             
	 _y.at(i)    += pos.y()*weight;             
	 _z.at(i)    += pos.z()*weight;             

	 _x2.at(i)   += pos.x()*pos.x()*weight;             
	 _y2.at(i)   += pos.y()*pos.y()*weight;             
	 _z2.at(i)   += pos.z()*pos.z()*weight;             

	 _xy.at(i)   += pos.x()*pos.y()*weight;             
	 _xz.at(i)   += pos.x()*pos.z()*weight;             
	 _yz.at(i)   += pos.y()*pos.z()*weight;             

	 _w.at(i)    += weight;             
	 _w2.at(i)   += weight*weight;             
    }


    void CaloShowerStepMCUtil::reset(int i)
    {
	_n[i]    = 0;
	_edep[i] = 0;
	_time[i] = 0;
	_x[i]    = _y[i]  = _z[i]  = 0;
	_x2[i]   = _y2[i] = _z2[i] = 0;
	_xy[i]   = _yz[i] = _xz[i] = 0;
	_w[i]    = _w2[i] = 0;	
    }

    
    CLHEP::Hep3Vector& CaloShowerStepMCUtil::pos(int i)
    {        
	_pos[0] = _x[i]/_w[i];
        _pos[1] = _y[i]/_w[i];
        _pos[2] = _z[i]/_w[i];
      
        return _pos;
    }

    CLHEP::HepSymMatrix& CaloShowerStepMCUtil::covPos(int i)
    {
        
        double norm = _w[i]/(_w[i]*_w[i]-_w2[i]);
      
        _cov[0][0] = norm*(_x2[i] - _x[i]*_x[i]/_w[i]);
        _cov[1][1] = norm*(_y2[i] - _y[i]*_y[i]/_w[i]);
        _cov[2][2] = norm*(_z2[i] - _z[i]*_z[i]/_w[i]);
        _cov[0][1] = norm*(_xy[i] - _x[i]*_y[i]/_w[i]);
        _cov[0][2] = norm*(_xz[i] - _x[i]*_z[i]/_w[i]);
        _cov[1][2] = norm*(_yz[i] - _y[i]*_z[i]/_w[i]);

        return _cov;
    }


    void CaloShowerStepMCUtil::printBucket(int i)
    {
         std::cout<<"Entries= "<<_n[i]<<" Energy = "<<_edep[i]<<" Time = "<<_time[i]/_w[i]<<" pos=("<<_x[i]/_w[i]<<","<<_y[i]/_w[i]<<","<<_z[i]/_w[i]<<")"<<std::endl;
    }

}



