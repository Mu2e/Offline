

#include "TrkReco/inc/TrkQualHelper.hh"

void mu2e::TrkQualHelper::fillTrkQual(const KalRep* krep, TrkQual& trkqual) {
    
  if (checkStatus(krep->fitStatus().success())) {
    std::vector<TrkStrawHitSeed> hits;
    TrkUtilities::fillHitSeeds(krep, hits);
    std::vector<TrkStraw> straws;
    TrkUtilities::fillStraws(krep, straws);
    fillHitCountVariables(hits, straws, trkqual);

    KalSegment kseg;
    findBestKSeg(krep, kseg);
    fillTrkVariables(TrkUtilities::chisqConsistency(krep), krep->t0(), kseg, trkqual);

    trkqual.setMVAValue(_trkqualmva->evalMVA(trkqual.values()));
    trkqual.setMVAStatus(TrkQual::calculated);
  }
}
  
void mu2e::TrkQualHelper::fillTrkQual(const KalSeed& kseed, TrkQual& trkqual) {
  if (checkStatus(kseed.status())) {
    fillHitCountVariables(kseed.hits(), kseed.straws(), trkqual);

    KalSegment kseg;
    bool foundBestKSeg = findBestKSeg(kseed, kseg);

    fillTrkVariables(kseed.fitConsistency(), kseed.t0(), kseg, trkqual);

    trkqual.setMVAValue(_trkqualmva->evalMVA(trkqual.values()));

    if (foundBestKSeg) {
      trkqual.setMVAStatus(TrkQual::calculated);
    }
    else {
      trkqual.setMVAStatus(TrkQual::filled);
    }
    //      trkqual.setMVAStatus(TrkQual::failed);
  }
}

bool mu2e::TrkQualHelper::checkStatus(const int& status) {
  if (status > 0) {
    return true;
  }
  else {
    return false;
  }
}

void mu2e::TrkQualHelper::findBestKSeg(const KalRep* krep, KalSegment& kseg) {
  double zpos = 0;
  double fltlen = krep->pieceTraj().zFlight(zpos);
  // sample the momentum at this flight.  This belongs in a separate utility FIXME
  BbrVectorErr momerr = krep->momentumErr(fltlen);
  // sample the helix
  double locflt(0.0);
  const HelixTraj* htraj = dynamic_cast<const HelixTraj*>(krep->localTrajectory(fltlen,locflt));
  // fill the segment
  TrkUtilities::fillSegment(*htraj, momerr, kseg);
}

bool mu2e::TrkQualHelper::findBestKSeg(const KalSeed& kseed, KalSegment& kseg) {
  std::vector<KalSegment> const& ksegs = kseed.segments();
  auto bestkseg = ksegs.begin();
  for(auto ikseg = ksegs.begin(); ikseg != ksegs.end(); ++ikseg){
    HelixVal const& hel = ikseg->helix();
    // check for a segment whose range includes z=0.  There should be a better way of doing this, FIXME
    double sind = hel.tanDip()/sqrt(1.0+hel.tanDip()*hel.tanDip());
    if(hel.z0()+sind*ikseg->fmin() < 0.0 && hel.z0()+sind*ikseg->fmax() > 0.0){
      bestkseg = ikseg;
      break;
    }
  }
  kseg = *bestkseg;
  if (bestkseg != ksegs.end()) {
    return true;
  }
  else {
    return false;
  }      
}

bool mu2e::TrkQualHelper::checkStatus(const TrkFitFlag& status) {
  static TrkFitFlag goodfit(TrkFitFlag::kalmanOK);
  if (status.hasAllProperties(goodfit)) {
    return true;
  }
  else {
    return false;
  }
}

void mu2e::TrkQualHelper::fillHitCountVariables(const std::vector<TrkStrawHitSeed>& hits, const std::vector<TrkStraw>& straws, TrkQual& trkqual) {
    
  unsigned nhits = 0; unsigned nactive = 0; unsigned ndactive = 0; unsigned nnullambig = 0;
  TrkUtilities::countHits(hits, nhits, nactive, ndactive, nnullambig);

  trkqual[TrkQual::nactive] = nactive;
  trkqual[TrkQual::factive] = (double)nactive / nhits;
  trkqual[TrkQual::fdouble] = (double)ndactive / nactive;
  trkqual[TrkQual::fnullambig] = (double)nnullambig / nactive;
  trkqual[TrkQual::fstraws] = (double)straws.size() / nactive;
}

void mu2e::TrkQualHelper::fillTrkVariables(const double& fit_con, const TrkT0& trkt0, const KalSegment& bestkseg, TrkQual& trkqual) {

  if (fit_con > FLT_MIN) {
    trkqual[TrkQual::log10fitcon] = log10(fit_con);
  }
  else {
    trkqual[TrkQual::log10fitcon] = -50.0;
  }
  trkqual[TrkQual::t0err] = trkt0.t0Err();
  trkqual[TrkQual::momerr] = bestkseg.momerr();
  trkqual[TrkQual::d0] = bestkseg.helix().d0();
  trkqual[TrkQual::rmax] = bestkseg.helix().d0() + 2.0/bestkseg.helix().omega();
}
