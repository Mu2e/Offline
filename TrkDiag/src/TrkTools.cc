//
// Namespace for collecting tools used in TrkDiag tree filling
// Original author: A. Edmonds (November 2018)
//
#include "TrkDiag/inc/TrkTools.hh"
#include "RecoDataProducts/inc/TrkStrawHitSeed.hh"

#include "TrackerGeom/inc/Tracker.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"

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

	const auto& jhit = ihit+1;
	const auto& hhit = ihit-1;
	if( (jhit != hits.end() &&
	     jhit->flag().hasAllProperties(active) &&
	     jhit->strawId().getPlane() == ihit->strawId().getPlane() &&
	     jhit->strawId().getPanel() == ihit->strawId().getPanel() ) ||
	    (hhit >= hits.begin() &&
	     hhit->flag().hasAllProperties(active) &&
	     hhit->strawId().getPlane() == ihit->strawId().getPlane() &&
	     hhit->strawId().getPanel() == ihit->strawId().getPanel() )
	    ) {
	  ++ndouble;
	  if (ihit->flag().hasAllProperties(StrawHitFlag::active)) {
	    ++ndactive;
	  }
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
      //    trkinfo._status = kseed.status(); //TODO
      trkinfo._pdg = kseed.particle().particleType();
      trkinfo._t0 = kseed.t0().t0();
      trkinfo._t0err = kseed.t0().t0Err();
      //    trkinfo._ndof = krep->nDof(); // TODO

      unsigned int nhits(-1), nactive(-1), ndouble(-1), ndactive(-1), nnullambig(-1);
      countHits(kseed.hits(), nhits, nactive, ndouble, ndactive, nnullambig);
      trkinfo._nhits = nhits;
      trkinfo._nactive = nactive;
      trkinfo._ndouble = ndouble;
      trkinfo._ndactive = ndactive;
      trkinfo._nnullambig = nnullambig;

      trkinfo._chisq = kseed.chisquared();
      trkinfo._fitcon = kseed.fitConsistency();

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

      //    trkinfo._nbend = nbend; // TODO
      const KalSegment& kseg = *(kseed.segments().begin()); // is this the correct segment to get? TODO
      trkinfo._ent._fitmom = kseg.mom();
      trkinfo._ent._fitmomerr = kseg.momerr();
      trkinfo._ent._fltlen = kseg.fmin(); //TODO
      trkinfo._ent._fitpar = kseg.helix();
      //      trkinfo._ent._fitparerr = kseg.covar(); //TODO
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
	// for now approximate the local bfield direction as the z axis FIXME!!
	tshinfo._bdot = tdir.z(); // TODO
	*/
	tshinfo._t0 = ihit->t0().t0();
	tshinfo._t0err = ihit->t0().t0Err(); //	was: tshinfo._t0err = ihit->t0Err()/ihit->driftVelocity();
	/*
	// include signal propagation time correction
	tshinfo._ht = ihit->time()-ihit->signalTime(); // TODO
	*/
	tshinfo._ambig = ihit->ambig();
	/*
	if(ihit->hasResidual())
	  tshinfo._doca = ihit->poca().doca(); // TODO
	else
	  tshinfo._doca = -100.0;
	tshinfo._exerr = ihit->driftVelocity()*ihit->temperature(); // TODO
	tshinfo._penerr = ihit->penaltyErr(); // TODO
	*/
	tshinfo._edep = ihit->energyDep();
	tshinfo._wdist = ihit->wireDist();
	tshinfo._werr = ihit->wireRes();	
	tshinfo._driftend = ihit->driftEnd();
	tshinfo._tdrift = ihit->driftTime(); // TODO


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
      // loop over sites, pick out the materials
      
      for(const auto& i_straw : kseed.straws()) { // TODO this isn't the correct thing to loop over
	TrkStrawMatInfo tminfo;

	tminfo._plane = i_straw.straw().getPlane();
	tminfo._panel = i_straw.straw().getPanel();
	tminfo._layer = i_straw.straw().getLayer();
	tminfo._straw = i_straw.straw().getStraw();

	tminfo._active = i_straw.active(); // TODO
	tminfo._dp = i_straw.pfrac();
	tminfo._radlen = i_straw.radLen();
	    /*
	    tminfo._sigMS = kmat->deflectRMS(); // TODO
	    // DetIntersection info
	    const DetIntersection& dinter = kmat->detIntersection();
	    tminfo._thit = (dinter.thit != 0); // TODO
	    tminfo._thita = (dinter.thit != 0 && dinter.thit->isActive()); // TODO
	    */
	tminfo._doca = i_straw.doca();
	tminfo._tlen = i_straw.trkLen();

	tminfos.push_back(tminfo);
      }
    }

    void fillCaloHitInfo(const KalSeed& kseed, TrkCaloHitInfo& tchinfo ) {

      if (kseed.hasCaloCluster()) {
	const TrkCaloHitSeed& tch = kseed.caloHit();

	tchinfo._active = tch.flag().hasAllProperties(StrawHitFlag::active);
	tchinfo._did = tch.caloCluster()->diskId();
	tchinfo._trklen = tch.trkLen();
	tchinfo._clen = tch.hitLen();
	//	  HepPoint hpos = tch->hitTraj()->position(tch->hitLen());
	tchinfo._poca = tch.caloCluster()->cog3Vector();//XYZVec(hpos.x(),hpos.y(),hpos.z()); //TODO
	if(tch.flag().hasAllProperties(StrawHitFlag::doca)) {
	  tchinfo._doca = tch.clusterAxisDOCA();
	}
	else {
	  tchinfo._doca = -100.0;
	}
	tchinfo._t0 = tch.t0().t0();
	tchinfo._t0err = tch.t0().t0Err();
	tchinfo._ct = tch.caloCluster()->time();
	tchinfo._edep = tch.caloCluster()->energyDep();
	//	  Hep3Vector tdir = tch->trkTraj()->direction(tchinfo._trklen);
	//tchinfo._cdot = tdir.dot(Hep3Vector(0.0,0.0,1.0)); // TODO

      }
    }
  }
}
