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
#include "ProditionsService/inc/ProditionsHandle.hh"
#include "Mu2eUtilities/inc/MVATools.hh"
#include "AnalysisConditions/inc/TrkQualCatalog.hh"
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
    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::Atom<art::InputTag> kalSeedTag{Name("KalSeedCollection"), Comment("Input tag for KalSeedCollection")};
      fhicl::Atom<std::string> trainName{Name("TrainingName"), Comment("Name of the training (e.g. TrkQual)")};
      fhicl::Atom<bool> printMVA{Name("PrintMVA"), Comment("Print the MVA used"), false};
    };

    using Parameters = art::EDProducer::Table<Config>;
    TrackQuality(const Parameters& conf);

  private:
    void produce(art::Event& event) override;
    void initializeMVA(std::string xmlfilename);

    art::InputTag _kalSeedTag;
    std::string _trainName;
    bool _printMVA;

    mu2e::ProditionsHandle<mu2e::TrkQualCatalog> _trkQualCatalogH;

    InfoStructHelper _infoStructHelper;
  };

  TrackQuality::TrackQuality(const Parameters& conf) :
    art::EDProducer{conf},
    _kalSeedTag(conf().kalSeedTag()), 
    _trainName(conf().trainName()),
    _printMVA(conf().printMVA())
  {
    produces<TrkQualCollection>();
    produces<RecoQualCollection>();
  }

  void TrackQuality::produce(art::Event& event ) {
    // create output
    unique_ptr<TrkQualCollection> tqcol(new TrkQualCollection());
    unique_ptr<RecoQualCollection> rqcol(new RecoQualCollection());

    // get the KalSeeds
    art::Handle<KalSeedCollection> kalSeedHandle;
    event.getByLabel(_kalSeedTag, kalSeedHandle);
    const auto& kalSeeds = *kalSeedHandle;

    TrkQualCatalog const& trkQualCatalog = _trkQualCatalogH.get(event.id());
    TrkQualEntry const& trkQualEntry = trkQualCatalog.find(_trainName);

    if(_printMVA) {
      trkQualEntry._mvaTool->showMVA();
    }

    // Go through the tracks and calculate their track qualities
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
	  trkqual.setMVAValue(trkQualEntry._mvaTool->evalMVA(trkqual.values(), trkQualEntry._mvaMask));

	}
	else {
	  trkqual.setMVAStatus(MVAStatus::filled);
	}
      }
      tqcol->push_back(trkqual);

      // Get the efficiency cut that this track passes
      Float_t passCalib = 0.0; // everything will pass a 100% efficient cut
      if (trkQualEntry._calibrated) {
	//	std::cout << _trainName << " = " << trkqual.MVAValue() << std::endl;
	for (const auto& i_pair : trkQualEntry._effCalib) {
	  //	  std::cout << i_pair.first << ", " << i_pair.second << ": ";
	  if (trkqual.MVAValue() >= i_pair.second) {
	    //	    std::cout << "PASSES" << std::endl;
	    passCalib = i_pair.first;
	  }
	  else {
	    //	    std::cout << "FAILS" << std::endl;
	    break;
	  }
	}
      }
      rqcol->push_back(RecoQual(trkqual.status(),trkqual.MVAValue(), passCalib));
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
