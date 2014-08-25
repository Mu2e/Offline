// utilities for the Module to perform BaBar Kalman fit
//
// $Id: TrkPatRecUtils.hh,v 1.2 2014/08/25 12:08:29 tassiell Exp $
// $Author: tassiell $
// $Date: 2014/08/25 12:08:29 $
//

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
// BaBar
#include "TrkPatRec/inc/TrkHitFilter.hh"
#include "TrkBase/TrkPoca.hh"
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
              tmpseed._fullTrkSeed._covMtrx[i][j]=seeddef.convMatr()[i+1][j+1];
            }
          }

          for (std::vector<hitIndex>::const_iterator ihit=seeddef.strawHitIndices().begin(); ihit!=seeddef.strawHitIndices().end(); ++ihit) {
                  tmpseed._fullTrkSeed._selectedTrackerHitsIdx.push_back( mu2e::HitIndex( ihit->_index, ihit->_ambig) );
                  tmpseed._selectedTrackerHits.push_back( StrawHitPtr (strawhitsH,ihit->_index) );
          }
  }

  void filterOutliers(TrkDef& mytrk,Trajectory const& traj,double maxdoca, int diag=0, vector<TrkHitFilter> *thfvec=0, KalDiag *kdiag=0);

}
