// utilities for the Module to perform BaBar Kalman fit
//
// $Id: TrkPatRecUtils.hh,v 1.1 2014/08/22 16:10:41 tassiell Exp $
// $Author: tassiell $
// $Date: 2014/08/22 16:10:41 $
//

// framework
#include "KalmanTests/inc/TrkDef.hh"
//#include "art/Framework/Principal/Handle.h"
#include "GeometryService/inc/GeomHandle.hh"
#include "art/Framework/Principal/Group.h"
// conditions
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TTrackerGeom/inc/TTracker.hh"
//#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
//Mu2e
#include "KalmanTests/inc/KalFitMC.hh"
#include "TrkPatRec/inc/TrkPatRec.hh"
// data
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StereoHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
//#include "RecoDataProducts/inc/TrackerHitTimeCluster.hh"
#include "RecoDataProducts/inc/TrackerHitTimeClusterCollection.hh"
#include "RecoDataProducts/inc/TrackSeed.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruth.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
// BaBar
#include "TrkPatRec/inc/TrkHitFilter.hh"
#include "TrkBase/TrkPoca.hh"
// root
#include "TH1F.h"
#include "TMVA/Reader.h"
//#include "TSpectrum.h"
// boost
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/accumulators/statistics/mean.hpp>
// C++
#include <vector>
#include <string>
#include <iostream>
#include "Rtypes.h"

using namespace std; 
using namespace boost::accumulators;

namespace mu2e 
{
  inline void cleanTimePeak(TrkTimePeak& tpeak,
                     const StrawHitCollection* shcol, const StrawHitPositionCollection* shpcol,
                     TimePeakMVA &pmva, TMVA::Reader *peakMVA, std::string &PMVAType,
                     double &minpeakmva, unsigned &minnhits,
                     double &maxpeakdt, double &maxpeakdphi) {
    static const double twopi(2*M_PI);
    // iteratively filter outliers
    double worstmva(100.0);
    double pphi(1000.0);
    double ptime(-1000.0);
    do{
  // first, compute the average phi and range.  Take care for wrapping
      accumulator_set<double, stats<tag::mean > > facc;
      accumulator_set<double, stats<tag::mean > > tacc;
      for(size_t ips=0;ips<tpeak._trkptrs.size();++ips){
        unsigned ish = tpeak._trkptrs[ips]._index;
        double time = shcol->at(ish).time();
        tacc(time);
        double phi = shpcol->at(ish).pos().phi();
        if(extract_result<tag::count>(facc) > 0){
          double dphi = phi - extract_result<tag::mean>(facc);
          if(dphi > M_PI){
            phi -= twopi;
          } else if(dphi < -M_PI){
            phi += twopi;
          }
        }
        facc(phi);
      }
      pphi = extract_result<tag::mean>(facc);
      ptime = extract_result<tag::mean>(tacc);
// find the least signal-like hit
      size_t iworst =0;
      worstmva = 100.0;
      for(size_t ips=0;ips<tpeak._trkptrs.size();++ips){
        unsigned ish = tpeak._trkptrs[ips]._index;
        double dt = shcol->at(ish).time() - ptime;
        double rho = shpcol->at(ish).pos().perp();
        double phi = shpcol->at(ish).pos().phi();
        double dphi = phi - pphi;
        if(dphi > M_PI){
          dphi -= twopi;
        } else if(dphi < -M_PI){
          dphi += twopi;
        }
        // compute MVA
        pmva._dt = dt;
        pmva._dphi = dphi;
        pmva._rho = rho;
        double mvaout = peakMVA->EvaluateMVA(PMVAType);
        if(mvaout < worstmva){
          worstmva = mvaout;
          iworst = ips;
        }
      }
      // remove the worst hit
      if(worstmva < minpeakmva){
        std::swap(tpeak._trkptrs[iworst],tpeak._trkptrs.back());
        tpeak._trkptrs.pop_back();
      }
    } while(tpeak._trkptrs.size() >= minnhits && worstmva < minpeakmva);
// final pass: hard cut on dt and dphi
    vector<size_t> toremove;
    accumulator_set<double, stats<tag::mean > > facc;
    accumulator_set<double, stats<tag::mean > > tacc;
    for(size_t ips=0;ips<tpeak._trkptrs.size();++ips){
      unsigned ish = tpeak._trkptrs[ips]._index;
      double dt = shcol->at(ish).time() - ptime;
      double phi = shpcol->at(ish).pos().phi();
      double dphi = phi - pphi;
      if(dphi > M_PI){
        dphi -= twopi;
        phi -= twopi;
      } else if(dphi < -M_PI){
        dphi += twopi;
        phi += twopi;
      }
      if(fabs(dt) < maxpeakdt && fabs(dphi) < maxpeakdphi){
        tacc(shcol->at(ish).time());
        facc(phi);
      } else {
        toremove.push_back(ips);
      }
    }
    // actually remove the hits; must start from the back
    std::sort(toremove.begin(),toremove.end(),std::greater<size_t>());
    for(auto irm=toremove.begin();irm!=toremove.end();++irm){
      std::swap(tpeak._trkptrs[*irm],tpeak._trkptrs.back());
      tpeak._trkptrs.pop_back();
    }
    // update peak properties
    pphi = extract_result<tag::mean>(facc);
    ptime = extract_result<tag::mean>(tacc);
    tpeak._tpeak = ptime;
    tpeak._phi = pphi;

  }

//  inline void initializeReaders(TMVA::Reader *peakMVA, TimePeakMVA &pmva, std::string &PMVAType, std::string &PMVAWeights) {
//    peakMVA = new TMVA::Reader();
//    peakMVA->AddVariable("_dt",&pmva._dt);
//    peakMVA->AddVariable("_dphi",&pmva._dphi);
//    peakMVA->AddVariable("_rho",&pmva._rho);
//    peakMVA->BookMVA(PMVAType,PMVAWeights);
//  }

  inline void findTimePeaks( TH1F const& tspect,std::vector<Float_t>& xpeak,std::vector<Float_t>& ypeak, double &onedthresh) {
    int ibin(1); // root starts bin numbers at 1
    int nbins = tspect.GetNbinsX();
    while(ibin < nbins+1){
      double y = tspect.GetBinContent(ibin);
      if(y >= onedthresh){
        double xmax = tspect.GetBinCenter(ibin);
        double ymax = y;
        double yprev = y;
        // find contiguous bins above threshold
        int jbin = ibin+1;
        bool descending(false);
        while(jbin < nbins+1 && tspect.GetBinContent(jbin) >= onedthresh){
          y =tspect.GetBinContent(jbin);
          descending |= yprev-y > sqrt(yprev);
          // don't follow next maximum
          if(descending && y-yprev > sqrt(y)){
            break;
          } else {
            if(y > ymax){
              ymax = y;
              xmax = tspect.GetBinCenter(jbin);
            }
            yprev = y;
            ibin = jbin;
            ++jbin;
          }
        }
        xpeak.push_back(xmax);
        ypeak.push_back(ymax);
      }
      ++ibin;
    }
  }

  inline void findTimePeaks(const StrawHitCollection* shcol, const StrawHitPositionCollection* shpcol, StrawHitFlagCollection* flags,
                  StrawHitFlag &tsel, StrawHitFlag &hsel, StrawHitFlag &tbkg, StrawHitFlag &hbkg,
                  std::vector<TrkTimePeak> &tpeaks, unsigned &maxnpeak,
                  unsigned &nbins, double &tmin, double &tmax,
                  double &onedthresh, double &ymin, double &maxdt, unsigned &minnhits,
                  bool &cleanpeaks, TimePeakMVA &pmva, TMVA::Reader *peakMVA, std::string &PMVAType, 
                  double &minpeakmva, double &maxpeakdt, double &maxpeakdphi) {
    TH1F timespec("timespec","time spectrum",nbins,tmin,tmax);
    // loop over straws hits and fill time spectrum plot for tight hits
    unsigned nstrs = shcol->size();
    for(unsigned istr=0; istr<nstrs;++istr){
      if(flags->at(istr).hasAllProperties(tsel) && !flags->at(istr).hasAnyProperty(tbkg)){
	double time = shcol->at(istr).time();
	timespec.Fill(time);
      }
    }
    std::vector<Float_t> xpeaks,ypeaks;
    findTimePeaks(timespec,xpeaks,ypeaks,onedthresh);
    unsigned np = xpeaks.size();
    // Loop over peaks, looking only at those with a minimum peak value
    for (unsigned ip=0; ip<np; ++ip) {
      Float_t xp = xpeaks[ip];
      Float_t yp = ypeaks[ip];
      TrkTimePeak tpeak(xp,yp);
      if(yp > ymin){
    // associate all hits in the time window with this peak
	for(size_t istr=0; istr<nstrs;++istr){
	  if(flags->at(istr).hasAllProperties(hsel) && !flags->at(istr).hasAnyProperty(hbkg)){
	    double time = shcol->at(istr).time();
	    if(fabs(time-xp) < maxdt){
	      tpeak._trkptrs.push_back(istr);
	    }
	  }
	}
	// if requested, clean the peaks
	if(cleanpeaks)cleanTimePeak(tpeak, shcol, shpcol, pmva, peakMVA, PMVAType, minpeakmva, minnhits, maxpeakdt, maxpeakdphi);
	if(tpeak._trkptrs.size() >= minnhits)tpeaks.push_back(tpeak);
      }
    }
    // sort the peaks so that the largest comes first
    sort(tpeaks.begin(),tpeaks.end(),greater<TrkTimePeak>());
//    // if requested, fill diagnostics
//    if(_diag>1 && _kfitmc.mcData()._mcsteps != 0){
//      for(size_t ip=0;ip<_tpeaks.size();++ip){
//	TrkTimePeak const& tp = _tpeaks[ip];
//	fillPeakDiag(ip,tp);
//      }
//    }
                    
  }

  inline void createTimePeak(const StrawHitCollection* shcol, StrawHitFlagCollection* flags,
                  StrawHitFlag &tsel, StrawHitFlag &hsel, StrawHitFlag &tbkg, StrawHitFlag &hbkg,
                  std::vector<TrkTimePeak> &tpeaks,
                  unsigned &minnhits) {
    // find the median time
    accumulator_set<double, stats<tag::median(with_p_square_quantile) > > tacc;
    unsigned nstrs = shcol->size();
    for(unsigned istr=0; istr<nstrs;++istr){
      if(flags->at(istr).hasAllProperties(tsel) && !flags->at(istr).hasAnyProperty(tbkg)){
	double time = shcol->at(istr).time();
	tacc(time);
      }
    }
    unsigned np = boost::accumulators::extract::count(tacc);  
    if(np >= minnhits){
      double mtime  = median(tacc);
// create a time peak from the full subset of selected hits
      TrkTimePeak tpeak(mtime,(double)nstrs);
      for(unsigned istr=0; istr<nstrs;++istr){
	if(flags->at(istr).hasAllProperties(hsel) && !flags->at(istr).hasAnyProperty(hbkg)){
	  tpeak._trkptrs.push_back(istr);
	}
      }
      tpeaks.push_back(tpeak);
    }
  }

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

  inline void filterOutliers(TrkDef& mytrk,Trajectory const& traj,double maxdoca, int diag=0, vector<TrkHitFilter> *thfvec=0, KalFitMC *kfitmc=0){
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
        if(kfitmc->mcData()._mcsteps != 0){
          const vector<MCHitSum>& mcsum = kfitmc->mcHitSummary(ihit);
          thfilter._mcpdg = mcsum[0]._pdgid;
          thfilter._mcgen = mcsum[0]._gid;
          thfilter._mcproc = mcsum[0]._pid;
        }
        thfvec->push_back(thfilter);
      }
    }
    // update track
    mytrk.setIndices(goodhits);
  }


}
