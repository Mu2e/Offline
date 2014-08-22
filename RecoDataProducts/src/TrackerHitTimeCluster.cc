//
// out data of the time peak algorithm for pattern recognition
//
// $Id: TrackerHitTimeCluster.cc,v 1.1 2014/08/22 16:10:41 tassiell Exp $
// $Author: tassiell $
// $Date: 2014/08/22 16:10:41 $
//
// Original author G. Tassielli
//

// Mu2e includes
#include "RecoDataProducts/inc/TrackerHitTimeCluster.hh"

// Root includes.
#include "TF1.h"

namespace mu2e {

  // Print the information found in this hit.
  void TrackerHitTimeCluster::print( ostream& ost, bool doEndl ) const {

    ost << *this;
    if ( doEndl ){
      ost << std::endl;
    }

  }

  void TrackerHitTimeCluster::expectedT0(double &t0, double &errt0, int type) const {
          switch (type) {
                case 3:
                {
                        t0 = _meanTime;
                        errt0 = _sigma;
                        break;
                }
                case 2:
                {
                        TF1 calib("toCalib","pol4",-200,50);
                        calib.SetParameters(0.577443,-0.00597028,4.23125e-06,3.33197e-07,1.17041e-09);
                        t0 = _minHitTime - calib.GetX((_meanTime-_minHitTime)/_nominalWidth);
                        errt0 = 0.16*_nominalWidth;
                        break;
                }
                case 1:
                {
                        if ( ((_meanTime-_minHitTime)/_nominalWidth)<0.8 ) {
                                t0 = (_minHitTime+(_meanTime-1.1*_sigma))*0.5;
                                errt0 = 0.092*_sigma;
                        } else {
                                t0 = _meanTime-0.98*_sigma;
                                errt0 = 0.23*_sigma;
                        }
                        break;
                }
                default:
                {
                        if ( ((_meanTime-_minHitTime)/_nominalWidth)<0.8 ) {
                                t0 = _minHitTime;
                                errt0 = 0.1*_sigma;
                        } else {
                                t0 = _meanTime-0.98*_sigma;
                                errt0 = 0.23*_sigma;
                        }
                        break;
                }
        }
  }

} // namespace mu2e
