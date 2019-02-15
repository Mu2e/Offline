//
// Namespace for collecting tools used in TrkDiag tree filling
// Original author: A. Edmonds (November 2018)
//
#include "TrkDiag/inc/TrkTools.hh"
#include "RecoDataProducts/inc/TrkStrawHitSeed.hh"

#include "TrackerGeom/inc/Tracker.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include <cmath>

namespace mu2e {
  namespace TrkTools {
    void countHits(const std::vector<TrkStrawHitSeed>& hits, unsigned& nhits, unsigned& nactive, unsigned& ndouble, unsigned& ndactive, unsigned& nnullambig) {
      nhits = 0; nactive = 0; ndouble = 0; ndactive = 0; nnullambig = 0;
      static StrawHitFlag active(StrawHitFlag::active);
      //      for (const auto& i_hit : hits) {
      for (std::vector<TrkStrawHitSeed>::const_iterator ihit = hits.begin(); ihit != hits.end(); ++ihit) {
	++nhits;
	if (ihit->flag().hasAllProperties(active)) {
	  ++nactive;
	  if (ihit->ambig()==0) {
	    ++nnullambig;
	  }
	}

	auto jhit = ihit; jhit++;
	if(jhit != hits.end() && ihit->strawId().uniquePanel() ==
	   jhit->strawId().uniquePanel()){
	    ++ndouble;
	    if(ihit->flag().hasAllProperties(active)) { ++ndactive; }
	}
      }
    }

    void countStraws(const std::vector<TrkStraw>& straws, unsigned& nmat, unsigned& nmatactive, double& radlen) {
      nmat = 0; nmatactive = 0; radlen = 0.0;
      for (std::vector<TrkStraw>::const_iterator i_straw = straws.begin(); i_straw != straws.end(); ++i_straw) {
	++nmat;
	if (i_straw->active()) {
	  ++nmatactive;
	  radlen += i_straw->radLen();
	}
      }
    }    

    void fillHitCount(StrawHitFlagCollection const& shfC, HitCount& hitcount) {
      hitcount._nsh = shfC.size();
      for(const auto& shf : shfC) {
	if(shf.hasAllProperties(StrawHitFlag::energysel))++hitcount._nesel;
	if(shf.hasAllProperties(StrawHitFlag::radsel))++hitcount._nrsel;
	if(shf.hasAllProperties(StrawHitFlag::timesel))++hitcount._ntsel;
	if(shf.hasAllProperties(StrawHitFlag::bkg))++hitcount._nbkg;
	if(shf.hasAllProperties(StrawHitFlag::stereo))++hitcount._nster;
	if(shf.hasAllProperties(StrawHitFlag::tdiv))++hitcount._ntdiv;
	if(shf.hasAllProperties(StrawHitFlag::trksel))++hitcount._ntpk;
	if(shf.hasAllProperties(StrawHitFlag::elecxtalk))++hitcount._nxt;
      }
    }


    void fillTrkInfo(const KalSeed& kseed,TrkInfo& trkinfo) {
      if(kseed.status().hasAllProperties(TrkFitFlag::kalmanConverged))
	trkinfo._status = 1;
      else if(kseed.status().hasAllProperties(TrkFitFlag::kalmanOK))
	trkinfo._status = 2;
      else
	trkinfo._status = -1;
      trkinfo._pdg = kseed.particle().particleType();
      trkinfo._t0 = kseed.t0().t0();
      trkinfo._t0err = kseed.t0().t0Err();

      unsigned int nhits(-1), nactive(-1), ndouble(-1), ndactive(-1), nnullambig(-1);
      countHits(kseed.hits(), nhits, nactive, ndouble, ndactive, nnullambig);
      trkinfo._ndof = nactive -5;
      trkinfo._nhits = nhits;
      trkinfo._nactive = nactive;
      trkinfo._ndouble = ndouble;
      trkinfo._ndactive = ndactive;
      trkinfo._nnullambig = nnullambig;

      trkinfo._chisq = kseed.chisquared();
      trkinfo._fitcon = kseed.fitConsistency();
      trkinfo._nbend = kseed.nBend();

      for(std::vector<TrkStrawHitSeed>::const_iterator ihit=kseed.hits().begin(); ihit != kseed.hits().end(); ++ihit) {
	if(ihit->flag().hasAllProperties(StrawHitFlag::active)) {
	  trkinfo._firstflt = ihit->trkLen();
	  break;
	}
      }
      for(std::vector<TrkStrawHitSeed>::const_reverse_iterator ihit=kseed.hits().rbegin(); ihit != kseed.hits().rend(); ++ihit) {
	if(ihit->flag().hasAllProperties(StrawHitFlag::active)) {
	  trkinfo._lastflt = ihit->trkLen();
	  break;
	}
      }

      // Loop through the KalSegments
      double firstflt = 9999999;
      double lastflt = -9999999;
      for (const auto& kseg : kseed.segments()) {
	//	std::cout << "AE: min = " << kseg.fmin() << ", max = " << kseg.fmax() << std::endl;
	if (kseg.fmin() < firstflt) {
	  firstflt = kseg.fmin();
	}
	if (kseg.fmax() > lastflt) {
	  lastflt = kseg.fmax();
	}
      }
      trkinfo._startvalid = firstflt;//krep->startValidRange();  // TODO
      trkinfo._endvalid = lastflt;//krep->endValidRange();  // TODO


      unsigned int nmat(-1), nmatactive(-1);
      double radlen(-1);
      countStraws(kseed.straws(), nmat, nmatactive, radlen);
      trkinfo._nmat = nmat;
      trkinfo._nmatactive = nmatactive;
      trkinfo._radlen = radlen; // TODO

      const KalSegment& kseg = *(kseed.segments().begin()); // is this the correct segment to get? TODO
      trkinfo._ent._fitmom = kseg.mom();
      trkinfo._ent._fitmomerr = kseg.momerr();
      trkinfo._ent._fltlen = kseg.fmin(); //TODO
      trkinfo._ent._fitpar = kseg.helix();
      CLHEP::HepSymMatrix pcov;
      kseg.covar().symMatrix(pcov);
      trkinfo._ent._fitparerr = helixpar(pcov);
    }

    void fillHitInfo(const KalSeed& kseed, std::vector<TrkStrawHitInfo>& tshinfos ) {
      tshinfos.clear();
      // loop over hits

      static StrawHitFlag active(StrawHitFlag::active);
      const Tracker& tracker = getTrackerOrThrow();
      for(std::vector<TrkStrawHitSeed>::const_iterator ihit=kseed.hits().begin(); ihit != kseed.hits().end(); ++ihit) {
	TrkStrawHitInfo tshinfo;

	tshinfo._active = ihit->flag().hasAllProperties(active);
	tshinfo._plane = ihit->strawId().plane();
	tshinfo._panel = ihit->strawId().panel();
	tshinfo._layer = ihit->strawId().layer();
	tshinfo._straw = ihit->strawId().straw();

	/*
	// kludge CLHEP problem
	HepPoint hpos = ihit->hitTraj()->position(ihit->hitLen());
	tshinfo._poca = XYZVec(hpos.x(),hpos.y(),hpos.z()); // TODO
	double resid,residerr;
	if(ihit->resid(resid,residerr,_uresid)){
	  tshinfo._resid = resid; // TODO
	  tshinfo._residerr = residerr; // TODO
	} else {
	  tshinfo._resid = tshinfo._residerr = -1000.;
	}
	*/
	tshinfo._rdrift = ihit->driftRad();
	tshinfo._rdrifterr = ihit->radialErr();
	
	double rstraw = tracker.getStraw(ihit->strawId()).getRadius();
	tshinfo._dx = std::sqrt(std::max(0.0,rstraw*rstraw-tshinfo._rdrift*tshinfo._rdrift));
	
	tshinfo._trklen = ihit->trkLen();
	tshinfo._hlen = ihit->hitLen();
	/*
	Hep3Vector tdir = ihit->trkTraj()->direction(tshinfo._trklen);
	tshinfo._wdot = tdir.dot(ihit->straw().getDirection()); // TODO
	*/
	tshinfo._t0 = ihit->t0().t0();
	tshinfo._t0err = ihit->t0().t0Err(); //	was: tshinfo._t0err = ihit->t0Err()/ihit->driftVelocity();
	tshinfo._ht = ihit->driftTime() + ihit->signalTime() + ihit->t0().t0();
	tshinfo._ambig = ihit->ambig();
	tshinfo._doca = ihit->wireDOCA();
	tshinfo._edep = ihit->energyDep();
	tshinfo._wdist = ihit->wireDist();
	tshinfo._werr = ihit->wireRes();	
	tshinfo._driftend = ihit->driftEnd();
	tshinfo._tdrift = ihit->driftTime(); 

	// count correlations with other TSH
	for(std::vector<TrkStrawHitSeed>::const_iterator jhit=kseed.hits().begin(); jhit != ihit; ++jhit) {
	  if(tshinfo._plane ==  jhit->strawId().plane() &&
	     tshinfo._panel == jhit->strawId().panel() ){
	    tshinfo._dhit = true;
	    if (jhit->flag().hasAllProperties(active)) {
	      tshinfo._dactive = true;
	      break;
	    }
	  }
	}

	tshinfos.push_back(tshinfo);
      }
    }

    void fillMatInfo(const KalSeed& kseed, std::vector<TrkStrawMatInfo>& tminfos ) {
      tminfos.clear();
     // loop over straws (material) 
      for(const auto& i_straw : kseed.straws()) { 
	TrkStrawMatInfo tminfo;

	tminfo._plane = i_straw.straw().getPlane();
	tminfo._panel = i_straw.straw().getPanel();
	tminfo._layer = i_straw.straw().getLayer();
	tminfo._straw = i_straw.straw().getStraw();

	tminfo._active = i_straw.active();
	tminfo._dp = i_straw.pfrac();
	tminfo._radlen = i_straw.radLen();
	    /*
	    tminfo._sigMS = kmat->deflectRMS(); // TODO
	    */
	tminfo._doca = i_straw.doca();
	tminfo._tlen = i_straw.trkLen();

	tminfos.push_back(tminfo);
      }
    }

    void fillCaloHitInfo(const KalSeed& kseed, 	Calorimeter const& calo, TrkCaloHitInfo& tchinfo) {
      if (kseed.hasCaloCluster()) {
	auto const& tch = kseed.caloHit();
	auto const& cc = tch.caloCluster();

	tchinfo._active = tch.flag().hasAllProperties(StrawHitFlag::active);
	tchinfo._did = cc->diskId();
	tchinfo._trklen = tch.trkLen();
	tchinfo._clen = tch.hitLen();
    
	if(tch.flag().hasAllProperties(StrawHitFlag::doca)) {
	  tchinfo._doca = tch.clusterAxisDOCA();
	}
	else {
	  tchinfo._doca = -100.0;
	}
	// add the propagation time offsetA
	tchinfo._t0 = tch.t0().t0();
	tchinfo._t0err = tch.t0().t0Err();
	tchinfo._ct = tch.time(); // time used to constrain T0 by this hit: includes the 'propagation time' offset
	tchinfo._cterr = tch.timeErr();
	tchinfo._edep = cc->energyDep();
	// transform cog to tracker coordinates; requires 2 steps.  This is at the front
	// of the disk
	XYZVec cpos = Geom::toXYZVec(calo.geomUtil().mu2eToTracker(calo.geomUtil().diskToMu2e(cc->diskId(),cc->cog3Vector())));
	// move to the front face and 
	// add the cluster length (relative to the front face).  crystal size should come from geom FIXME!
	cpos.SetZ(cpos.z() -200.0 + tch.hitLen());
	tchinfo._poca = cpos;
      }
    }

    void fillTrkQualInfo(const TrkQual& tqual, TrkQualInfo& trkqualInfo) {
      int n_trkqual_vars = TrkQual::n_vars;
      for (int i_trkqual_var = 0; i_trkqual_var < n_trkqual_vars; ++i_trkqual_var) {
	TrkQual::MVA_varindex i_index = TrkQual::MVA_varindex(i_trkqual_var);
	trkqualInfo._trkqualvars[i_trkqual_var] = (double) tqual[i_index];
      }
      trkqualInfo._trkqual = tqual.MVAOutput();
    }
  }
}
