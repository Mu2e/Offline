// utilities for the Module to perform BaBar Kalman fit
//
// $Id: TrkPatRecUtils.hh,v 1.3 2014/08/30 12:19:38 tassiell Exp $
// $Author: tassiell $
// $Date: 2014/08/30 12:19:38 $
//
// Original author G. Tassielli
//
#ifndef TrkPatRecUtils_HH
#define TrkPatRecUtils_HH

// framework
#include "KalmanTests/inc/TrkDef.hh"
//Mu2e
#include "KalmanTests/inc/KalDiag.hh"
#include "TrkPatRec/inc/TrkPatRec.hh"
// data
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/TrackerHitTimeClusterCollection.hh"
#include "RecoDataProducts/inc/TrackSeed.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "RecoDataProducts/inc/HelixVal.hh"
// BaBar
#include "BTrk/TrkBase/HelixTraj.hh"
#include "TrkPatRec/inc/TrkHitFilter.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
// C++
#include <vector>

namespace mu2e 
{

  inline void loadTimePeaks(std::vector<TrkTimePeak> &tpeaks, const TrackerHitTimeClusterCollection* tccol) {

    tpeaks.clear();
    for (size_t ipeak=0; ipeak<tccol->size(); ipeak++) {
      TrackerHitTimeCluster const&  tclust(tccol->at(ipeak));
      TrkTimePeak tpeak(tclust._meanTime,tclust._peakmax);
      for (std::vector<StrawHitPtr>::const_iterator thit=tclust._selectedTrackerHits.begin(); thit!=tclust._selectedTrackerHits.end(); ++thit) {
        tpeak._trkptrs.push_back(thit->key());
      }
      tpeaks.push_back(tpeak);
    }

    std::sort(tpeaks.begin(),tpeaks.end(),greater<TrkTimePeak>());
  }

  inline void fillTrackSeed(TrackSeed &tmpseed, TrkDef &seeddef,  art::Handle<TrackerHitTimeClusterCollection> &tclusthitH,  unsigned &ipeak, art::Handle<mu2e::StrawHitCollection> &strawhitsH) {

          tmpseed._relatedTimeCluster=TrackerHitTimeClusterPtr(tclusthitH,ipeak);
          tclusthitH.product()->at(ipeak).expectedT0(tmpseed._t0,tmpseed._errt0,3);
          tmpseed._fullTrkSeed._d0=seeddef.helix().d0();//_dpar[ParIndex::d0Index];//
          tmpseed._fullTrkSeed._phi0=seeddef.helix().phi0();
          tmpseed._fullTrkSeed._omega=seeddef.helix().omega();
          tmpseed._fullTrkSeed._z0=seeddef.helix().z0();
          tmpseed._fullTrkSeed._tanDip=seeddef.helix().tanDip();
          for(int i=0;i<5;i++) {
            for(int j=0;j<5;j++) {
              tmpseed._fullTrkSeed._covMtrx[i][j]=seeddef.helixCovMatr()[i][j];
            }
          }

          for (std::vector<hitIndex>::const_iterator ihit=seeddef.strawHitIndices().begin(); ihit!=seeddef.strawHitIndices().end(); ++ihit) {
                  tmpseed._fullTrkSeed._selectedTrackerHitsIdx.push_back( mu2e::HitIndex( ihit->_index, ihit->_ambig) );
                  tmpseed._selectedTrackerHits.push_back( StrawHitPtr (strawhitsH,ihit->_index) );
          }
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

