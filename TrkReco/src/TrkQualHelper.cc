

#include "TrkReco/inc/TrkQualHelper.hh"

 
void mu2e::TrkQualHelper::fillTrkQual(const KalSeed& kseed, TrkQual& trkqual, const double& zpos) {

  static TrkFitFlag goodfit(TrkFitFlag::kalmanOK);
  if (kseed.status().hasAllProperties(goodfit)) {

    // fill the hit count variables
    unsigned nhits = 0; unsigned nactive = 0; unsigned ndouble = 0; unsigned ndactive = 0; unsigned nnullambig = 0;
    TrkUtilities::countHits(kseed.hits(), nhits, nactive, ndouble, ndactive, nnullambig);
    trkqual[TrkQual::nactive] = nactive;
    trkqual[TrkQual::factive] = (double)nactive / nhits;
    trkqual[TrkQual::fdouble] = (double)ndactive / nactive;
    trkqual[TrkQual::fnullambig] = (double)nnullambig / nactive;
    trkqual[TrkQual::fstraws] = (double)kseed.straws().size() / nactive;

    // fill fit consistency and t0 variables
    if (kseed.fitConsistency() > FLT_MIN) {
      trkqual[TrkQual::log10fitcon] = log10(kseed.fitConsistency());
    }
    else {
      trkqual[TrkQual::log10fitcon] = -50.0;
    }
    trkqual[TrkQual::t0err] = kseed.t0().t0Err();

    // find the best KalSegment and fill variables relating to it
    KalSegment kseg;
    std::vector<KalSegment> const& ksegs = kseed.segments();
    auto bestkseg = ksegs.begin();
    for(auto ikseg = ksegs.begin(); ikseg != ksegs.end(); ++ikseg){
      HelixVal const& hel = ikseg->helix();
      // check for a segment whose range includes z=0.  There should be a better way of doing this, FIXME
      double sind = hel.tanDip()/sqrt(1.0+hel.tanDip()*hel.tanDip());
      if(hel.z0()+sind*ikseg->fmin() < zpos && hel.z0()+sind*ikseg->fmax() > zpos){
	bestkseg = ikseg;
	break;
      }
    }
    kseg = *bestkseg;
    if (bestkseg != ksegs.end()) {
      double charge = kseed.particle().charge();

      trkqual[TrkQual::momerr] = bestkseg->momerr();
      trkqual[TrkQual::d0] = -1*charge*bestkseg->helix().d0();
      trkqual[TrkQual::rmax] = -1*charge*(bestkseg->helix().d0() + 2.0/bestkseg->helix().omega());

      trkqual.setMVAStatus(TrkQual::calculated);
      trkqual.setMVAValue(_trkqualmva->evalMVA(trkqual.values()));
    }
    else {
      trkqual.setMVAStatus(TrkQual::filled);
    }
  }
}
