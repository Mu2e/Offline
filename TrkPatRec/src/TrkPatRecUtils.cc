// utilities for the Module to perform BaBar Kalman fit
//
// $Id: TrkPatRecUtils.cc,v 1.1 2014/08/25 12:08:29 tassiell Exp $
// $Author: tassiell $
// $Date: 2014/08/25 12:08:29 $
//

// framework
#include "art/Framework/Principal/Handle.h"
#include "GeometryService/inc/GeomHandle.hh"
// conditions
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TTrackerGeom/inc/TTracker.hh"
// data
#include "RecoDataProducts/inc/StrawDigiCollection.hh"
#include "MCDataProducts/inc/StrawDigiMCCollection.hh"
// Mu2e
#include "TrkPatRec/inc/TrkPatRecUtils.hh"

namespace mu2e 
{

  void filterOutliers(TrkDef& mytrk,Trajectory const& traj,double maxdoca, int diag, vector<TrkHitFilter> *thfvec, KalDiag *kdiag){
    //  Trajectory info
    Hep3Vector tdir;
    HepPoint tpos;
    traj.getInfo(0.0,tpos,tdir);
    // tracker and conditions
    const Tracker& tracker = getTrackerOrThrow();
    ConditionsHandle<TrackerCalibrations> tcal("ignored");
    const StrawHitCollection* hits = mytrk.strawHitCollection();
    const vector<hitIndex>& indices = mytrk.strawHitIndices();
    vector<hitIndex> goodhits;
    for(unsigned ihit=0;ihit<indices.size();++ihit){
      StrawHit const& sh = hits->at(indices[ihit]._index);
      Straw const& straw = tracker.getStraw(sh.strawIndex());
      CLHEP::Hep3Vector hpos = straw.getMidPoint();
      CLHEP::Hep3Vector hdir = straw.getDirection();
      // convert to HepPoint to satisfy antique BaBar interface: FIXME!!!
      HepPoint spt(hpos.x(),hpos.y(),hpos.z());
      TrkLineTraj htraj(spt,hdir,-20,20);
      // estimate flightlength along track.  This assumes a constant BField!!!
      double fltlen = (hpos.z()-tpos.z())/tdir.z();
      TrkPoca hitpoca(traj,fltlen,htraj,0.0);
      // flag hits with small residuals
      if(fabs(hitpoca.doca()) < maxdoca){
        goodhits.push_back(indices[ihit]);
      }
      // optional diagnostics
      if(diag > 0){
              // summarize the MC truth for this strawhit
              TrkHitFilter thfilter;
              HepPoint tpos =  traj.position(hitpoca.flt1());
              thfilter._pos = CLHEP::Hep3Vector(tpos.x(),tpos.y(),tpos.z());
              thfilter._doca = hitpoca.doca();
              if(kdiag->mcData()._mcdigis != 0){
                StrawDigiMC const& mcdigi = kdiag->mcData()._mcdigis->at(indices[ihit]._index);
                // use TDC channel 0 to define the MC match
                StrawDigi::TDCChannel itdc = StrawDigi::zero;
                if(!mcdigi.hasTDC(StrawDigi::one)) itdc = StrawDigi::one;
                art::Ptr<StepPointMC> const& spmcp = mcdigi.stepPointMC(itdc);
                art::Ptr<SimParticle> const& spp = spmcp->simParticle();
                thfilter._mcpdg = spp->pdgId();
                thfilter._mcproc = spp->creationCode();
                thfilter._mcgen = -1;
                if(spp->genParticle().isNonnull())
                  thfilter._mcgen = spp->genParticle()->generatorId().id();
              }
              thfvec->push_back(thfilter);
            }
    }
    // update track
    mytrk.setIndices(goodhits);
  }


}
