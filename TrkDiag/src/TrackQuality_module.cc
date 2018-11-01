//
// Create a TrkQual object
//
// Original author A. Edmonds
//

// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Utilities/make_tool.h"
// utilities
#include "Mu2eUtilities/inc/MVATools.hh"
#include "TrkReco/inc/TrkUtilities.hh"
// data
#include "RecoDataProducts/inc/KalSeed.hh"
#include "RecoDataProducts/inc/TrkQual.hh"
// C++
#include <iostream>
#include <fstream>
#include <string>
#include <functional>
#include <float.h>
#include <vector>
using namespace std;
using CLHEP::Hep3Vector;
using CLHEP::HepVector;

namespace mu2e
{

  class TrackQuality : public art::EDProducer
  {
  public:
    TrackQuality(fhicl::ParameterSet const&);

  private:
    void produce(art::Event& event) override;

    art::InputTag _kalSeedTag;

    MVATools* _trkqualmva;
  };

  TrackQuality::TrackQuality(fhicl::ParameterSet const& pset) :
    _kalSeedTag(pset.get<art::InputTag>("KalSeedCollection", "")),
    _trkqualmva(new MVATools(pset.get<fhicl::ParameterSet>("TrkQualMVA", fhicl::ParameterSet())))

  {
    produces<TrkQualCollection>();
    
    _trkqualmva->initMVA();
  }

  void TrackQuality::produce(art::Event& event ) {
    // create output
    unique_ptr<TrkQualCollection> tqcol(new TrkQualCollection());

    // get the KalSeeds
    art::Handle<KalSeedCollection> kalSeedHandle;
    event.getByLabel(_kalSeedTag, kalSeedHandle);
    const auto& kalSeeds = *kalSeedHandle;

    for (const auto& i_kalSeed : kalSeeds) {
      TrkQual trkqual;

      static TrkFitFlag goodfit(TrkFitFlag::kalmanOK);
      if (i_kalSeed.status().hasAllProperties(goodfit)) {

	// fill the hit count variables
	unsigned nhits = 0; unsigned nactive = 0; unsigned ndouble = 0; unsigned ndactive = 0; unsigned nnullambig = 0;
	TrkUtilities::countHits(i_kalSeed.hits(), nhits, nactive, ndouble, ndactive, nnullambig);
	trkqual[TrkQual::nactive] = nactive;
	trkqual[TrkQual::factive] = (double)nactive / nhits;
	trkqual[TrkQual::fdouble] = (double)ndactive / nactive;
	trkqual[TrkQual::fnullambig] = (double)nnullambig / nactive;
	trkqual[TrkQual::fstraws] = (double)i_kalSeed.straws().size() / nactive;

	// fill fit consistency and t0 variables
	if (i_kalSeed.fitConsistency() > FLT_MIN) {
	  trkqual[TrkQual::log10fitcon] = log10(i_kalSeed.fitConsistency());
	}
	else {
	  trkqual[TrkQual::log10fitcon] = -50.0;
	}
	trkqual[TrkQual::t0err] = i_kalSeed.t0().t0Err();
	
	// find the best KalSegment and fill variables relating to it
	KalSegment kseg;
	std::vector<KalSegment> const& ksegs = i_kalSeed.segments();
	auto bestkseg = ksegs.begin();
	double zpos = -1631.11;
	for(auto ikseg = ksegs.begin(); ikseg != ksegs.end(); ++ikseg){
	  HelixVal const& hel = ikseg->helix();
	  // check for a segment whose range includes zpos.  There should be a better way of doing this, FIXME
	  double sind = hel.tanDip()/sqrt(1.0+hel.tanDip()*hel.tanDip());
	  if(hel.z0()+sind*ikseg->fmin() < zpos && hel.z0()+sind*ikseg->fmax() > zpos){
	    bestkseg = ikseg;
	    break;
	  }
	}
	kseg = *bestkseg;
	if (bestkseg != ksegs.end()) {
	  double charge = i_kalSeed.particle().charge();
	  
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

      tqcol->push_back(trkqual);
    }

    // put the output products into the event
    event.put(move(tqcol));
  }  
}// mu2e

DEFINE_ART_MODULE(mu2e::TrackQuality);
