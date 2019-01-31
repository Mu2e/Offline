//
// Namespace for collecting tools used in TrkDiag tree filling
// Original author: A. Edmonds (November 2018)
//
#include "TrkDiag/inc/TrkTools.hh"
#include "RecoDataProducts/inc/TrkStrawHitSeed.hh"

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
	  /*	  if (ihit->nStrawHits()>=2) {
	    ++ndactive;
	  }
	  */
	  //	  std::cout << "AE: ihit->nStrawHits() = " << ihit->nStrawHits() << std::endl;
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
      //    trkinfo._radlen = krep->radiationFraction(); // TODO
      //trkinfo._firstflt = // TODO
      //trkinfo._lastflt = // TODOa
      //    trkinfo._startvalid = krep->startValidRange();  // TODO
      //    trkinfo._endvalid = krep->endValidRange();  // TODO
      //    trkinfo._nmat = nmat;  // TODO
      //    trkinfo._nmatactive = nmatactive; // TODO
      //    trkinfo._nbend = nbend; // TODO
      const KalSegment& kseg = *(kseed.segments().begin()); // is this the correct segment to get? TODO
      trkinfo._ent._fitmom = kseg.mom();
      trkinfo._ent._fitmomerr = kseg.momerr();
      //trkinfo.fltlen = ; //TODO
      trkinfo._ent._fitpar = kseg.helix();
      //    trkinfo._ent._fitparerr = kseg.covar(); //TODO
    }

    void fillHitInfo(const KalSeed& kseed, std::vector<TrkStrawHitInfo>& tshinfos ) {
      tshinfos.clear();
      // loop over hits

      static StrawHitFlag active(StrawHitFlag::active);
      for(std::vector<TrkStrawHitSeed>::const_iterator ihit=kseed.hits().begin(); ihit != kseed.hits().end(); ++ihit) {
	TrkStrawHitInfo tshinfo;

	tshinfo._active = ihit->flag().hasAllProperties(active);
	tshinfo._plane = ihit->strawId().plane();
	tshinfo._panel = ihit->strawId().panel();
	tshinfo._layer = ihit->strawId().layer();
	tshinfo._straw = ihit->strawId().straw();
	/*
	tshinfo._edep = ihit->comboHit().energyDep();
	// kludge CLHEP problem
	HepPoint hpos = ihit->hitTraj()->position(ihit->hitLen());
	tshinfo._poca = XYZVec(hpos.x(),hpos.y(),hpos.z());
	double resid,residerr;
	if(ihit->resid(resid,residerr,_uresid)){
	  tshinfo._resid = resid;
	  tshinfo._residerr = residerr;
	} else {
	  tshinfo._resid = tshinfo._residerr = -1000.;
	}
	tshinfo._rdrift = ihit->driftRadius();
	tshinfo._rdrifterr = ihit->driftRadiusErr();
	tshinfo._wdist = ihit->comboHit().wireDist();
	tshinfo._werr = ihit->comboHit().wireRes();
	double rstraw = ihit->straw().getRadius();
	tshinfo._dx = sqrt(max(0.0,rstraw*rstraw-tshinfo._rdrift*tshinfo._rdrift));
	tshinfo._trklen = ihit->fltLen();
	tshinfo._hlen = ihit->hitLen();
	Hep3Vector tdir = ihit->trkTraj()->direction(tshinfo._trklen);
	tshinfo._wdot = tdir.dot(ihit->straw().getDirection());
	// for now approximate the local bfield direction as the z axis FIXME!!
	tshinfo._bdot = tdir.z();
	tshinfo._t0 = ihit->hitT0()._t0;
	// include signal propagation time correction
	tshinfo._ht = ihit->time()-ihit->signalTime();
	//    tshinfo._dt = ihit->comboHit().dt();
	tshinfo._ambig = ihit->ambig();
	tshinfo._driftend = ihit->comboHit().driftEnd();
	tshinfo._tdrift = ihit->driftTime();
	//    tshinfo._totcal = ihit->comboHit().TOT(StrawEnd::cal);
	//    tshinfo._tothv = ihit->comboHit().TOT(StrawEnd::hv);
	if(ihit->hasResidual())
	  tshinfo._doca = ihit->poca().doca();
	else
	  tshinfo._doca = -100.0;
	tshinfo._exerr = ihit->driftVelocity()*ihit->temperature();
	tshinfo._penerr = ihit->penaltyErr();
	tshinfo._t0err = ihit->t0Err()/ihit->driftVelocity();
	// cannot count correlations with other hits in this function; set to false
	tshinfo._dhit = tshinfo._dactive = false;
	
	// count correlations with other TSH
	for(auto jhit=tshv.begin(); jhit != ihit; ++jhit){
	  const TrkStrawHit* otsh = *jhit;
	  if(otsh != 0){
	    if(tshinfo._plane ==  otsh->straw().id().getPlane() &&
	       tshinfo._panel == otsh->straw().id().getPanel() ){
	      tshinfo._dhit = true;
	      if(otsh->isActive()){
		tshinfo._dactive = true;
		break;
	      }
	    }
	  }
	}
	*/
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
	/*	if(isite->kalMaterial() != 0){
	  TrkStrawMatInfo tminfo;
	  const DetStrawElem* delem = dynamic_cast<const DetStrawElem*>(kmat->detElem());
	  if(delem != 0){
	    retval = true;
	    // KalMaterial info
	    tminfo._active = kmat->isActive();
	    tminfo._dp = kmat->momFraction();
	    tminfo._radlen = kmat->radiationFraction();
	    tminfo._sigMS = kmat->deflectRMS();
	    // DetIntersection info
	    const DetIntersection& dinter = kmat->detIntersection();
	    tminfo._thit = (dinter.thit != 0);
	    tminfo._thita = (dinter.thit != 0 && dinter.thit->isActive());
	    tminfo._doca = dinter.dist;
	    tminfo._tlen = dinter.pathlen;
	    // straw information
	    Straw const* straw = delem->straw();
	    
	    tminfos.push_back(tminfo);
	  }
	}
	*/
	tminfos.push_back(tminfo);
      }
    }
  }
}
