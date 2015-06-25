//
// Utilities functions to manage track seeds
//
// $Id: TrackSeedUtils.hh,v 1.2 2012/05/23 12:26:05 ignatov Exp $
// $Author: ignatov $
// $Date: 2012/05/23 12:26:05 $
//
// Original author G. Tassielli
//
#ifndef TrackSeedUtils_HH
#define TrackSeedUtils_HH

// Mu2e includes.
#include "RecoDataProducts/inc/HelixVal.hh"
#include "RecoDataProducts/inc/TrackSeed.hh"

// BaBar
#include "BTrk/BaBar/BaBar.hh"
//#include "KalmanTests/inc/TrkDef.hh"

//using namespace std;

namespace mu2e {

inline void HelixTraj2HelixVal (const HelixTraj & helIn, HelixVal &helOut) {
        helOut._d0     = helIn.d0();
        helOut._phi0   = helIn.phi0();
        helOut._omega  = helIn.omega();
        helOut._z0     = helIn.z0();
        helOut._tanDip = helIn.tanDip();
}

inline void HelixVal2HelixTraj (const HelixVal &helIn, HelixTraj &helOut) {
        //TrkExchangePar helParams( helIn._d0, helIn._phi0, helIn._omega, helIn._z0, helIn._tanDip );
        HepVector helParams(5);
        helParams(1) = helIn._d0;
        helParams(2) = helIn._phi0;
        helParams(3) = helIn._omega;
        helParams(4) = helIn._z0;
        helParams(5) = helIn._tanDip;
        HepSymMatrix conv(5,1);
        conv(1,1)=helIn._covMtrx[0][0]; conv(1,2)=helIn._covMtrx[0][1]; conv(1,3)=helIn._covMtrx[0][2]; conv(1,4)=helIn._covMtrx[0][3]; conv(1,5)=helIn._covMtrx[0][4];
        conv(2,1)=helIn._covMtrx[1][0]; conv(2,2)=helIn._covMtrx[1][1]; conv(2,3)=helIn._covMtrx[1][2]; conv(2,4)=helIn._covMtrx[1][3]; conv(2,5)=helIn._covMtrx[1][4];
        conv(3,1)=helIn._covMtrx[2][0]; conv(3,2)=helIn._covMtrx[2][1]; conv(3,3)=helIn._covMtrx[2][2]; conv(3,4)=helIn._covMtrx[2][3]; conv(3,5)=helIn._covMtrx[2][4];
        conv(4,1)=helIn._covMtrx[3][0]; conv(4,2)=helIn._covMtrx[3][1]; conv(4,3)=helIn._covMtrx[3][2]; conv(4,4)=helIn._covMtrx[3][3]; conv(4,5)=helIn._covMtrx[3][4];
        conv(5,1)=helIn._covMtrx[4][0]; conv(5,2)=helIn._covMtrx[4][1]; conv(5,3)=helIn._covMtrx[4][2]; conv(5,4)=helIn._covMtrx[4][3]; conv(5,5)=helIn._covMtrx[4][4];
        //helParams.setError(conv);
        HelixTraj tmpHelix(helParams,conv);
        helOut=tmpHelix;
}

}  // end namespace mu2e
#endif
