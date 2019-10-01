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
#include "art_root_io/TFileService.h"
#include "art/Utilities/make_tool.h"
// utilities
#include "Mu2eUtilities/inc/MVATools.hh"
#include "TrkDiag/inc/InfoStructHelper.hh"
#include "TrkDiag/inc/TrkInfo.hh"
// data
#include "RecoDataProducts/inc/KalSeed.hh"
#include "RecoDataProducts/inc/TrkQual.hh"
#include "RecoDataProducts/inc/RecoQual.hh"
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
    MVAMask _mvamask;

    InfoStructHelper _infoStructHelper;
  };

  TrackQuality::TrackQuality(fhicl::ParameterSet const& pset) :
    art::EDProducer{pset},
    _kalSeedTag(pset.get<art::InputTag>("KalSeedCollection", "")),
    _trkqualmva(new MVATools(pset.get<fhicl::ParameterSet>("TrkQualMVA", fhicl::ParameterSet())))

  {
    produces<TrkQualCollection>();
    produces<RecoQualCollection>();
    
    _trkqualmva->initMVA();

    // create the MVA mask in case we have removed variables
    const auto& labels = _trkqualmva->labels();
    _mvamask = 0;
    for (int i_var = 0; i_var < TrkQual::n_vars; ++i_var) {
      for (const auto& i_label : labels) {
	std::string i_varName = TrkQual::varName(static_cast<TrkQual::MVA_varindex>(i_var));
	if (i_label.find(i_varName) != std::string::npos) {
	  _mvamask ^= (1 << i_var);
	  break;
	}
      }
    }
    if(pset.get<bool>("PrintMVA",false))
      _trkqualmva->showMVA();
  }

  void TrackQuality::produce(art::Event& event ) {
    // create output
    unique_ptr<TrkQualCollection> tqcol(new TrkQualCollection());
    unique_ptr<RecoQualCollection> rqcol(new RecoQualCollection());

    // get the KalSeeds
    art::Handle<KalSeedCollection> kalSeedHandle;
    event.getByLabel(_kalSeedTag, kalSeedHandle);
    const auto& kalSeeds = *kalSeedHandle;

    for (const auto& i_kalSeed : kalSeeds) {
      TrkQual trkqual;

      static TrkFitFlag goodfit(TrkFitFlag::kalmanOK);
      if (i_kalSeed.status().hasAllProperties(goodfit)) {

	// fill the hit count variables
	TrkInfo tinfo;
	_infoStructHelper.fillTrkInfoHits(i_kalSeed, tinfo);
	_infoStructHelper.fillTrkInfoStraws(i_kalSeed, tinfo);
	trkqual[TrkQual::nactive] = tinfo._nactive;
	trkqual[TrkQual::factive] = (double) tinfo._nactive / tinfo._nhits;
	trkqual[TrkQual::fdouble] = (double) tinfo._ndactive / tinfo._nactive;
	trkqual[TrkQual::fnullambig] = (double) tinfo._nnullambig / tinfo._nactive;
	trkqual[TrkQual::fstraws] = (double)tinfo._nmatactive / tinfo._nactive;

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
	  
	  trkqual.setMVAStatus(MVAStatus::calculated);
	  trkqual.setMVAValue(_trkqualmva->evalMVA(trkqual.values(), _mvamask));

	}
	else {
	  trkqual.setMVAStatus(MVAStatus::filled);
	}
      }
      tqcol->push_back(trkqual);
      rqcol->push_back(RecoQual(trkqual.status(),trkqual.MVAValue()));
    }

    if ( (tqcol->size() != rqcol->size()) || (tqcol->size() != kalSeeds.size()) ) {
	throw cet::exception("TrackQuality") << "KalSeed, TrkQual and RecoQual sizes are inconsistent (" << kalSeeds.size() << ", " << tqcol->size() << ", " << rqcol->size() << " respectively)";
      }

    // put the output products into the event
    event.put(move(tqcol));
    event.put(move(rqcol));
  }  
}// mu2e

DEFINE_ART_MODULE(mu2e::TrackQuality);
