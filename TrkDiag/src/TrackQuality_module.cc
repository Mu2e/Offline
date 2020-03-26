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
#include "TrkDiag/inc/InfoStructHelper.hh"
#include "TrkDiag/inc/TrkInfo.hh"
#include "TrkDiag/inc/TrkQualCatalog.hh"
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

    //    mu2e::DbHandle<mu2e::TrkQualDb> _trkQualDb;
    mu2e::ProditionsHandle<mu2e::TrkQualCatalog> _trkQualCatalogH;
    //    std::string _currentXmlFile;
    //    MVATools* _trkqualmva;
    //    MVAMask _mvamask;
    std::map<Float_t, Float_t> _effCalib;

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
    for (const auto& i_trkQualEntry : trkQualCatalog.entries()) {
      if (i_trkQualEntry._trainName == _trainName) {

	if(_printMVA) {
	  i_trkQualEntry._mvaTool->showMVA();
	}


	// // get the XML filename for this TrkQual training
	// auto const& trkQualTable = _trkQualDb.get(event.id());
	// for (const auto& i_row : trkQualTable.rows()) {
	//   if (i_row.mvaname() == _trainName) { // check the training name
	// 	std::string xmlfilename = i_row.xmlfilename();
	// 	if (xmlfilename != _currentXmlFile) { // only reinitialize if the XML file has changed
	// 	  initializeMVA(xmlfilename);

	// 	  if (i_row.calibrated() == 1) {
	// 	    _effCalib = trkQualTable.getCalib(i_row.idx()); // get the calibration if it exists
	// 	  }
	// 	}
	// 	break;
	//   }
	// }
	// if(!i_trkQualEntry._mvaTool) {
	//   throw cet::exception("TrackQuality") << "i_trkQualEntry._mvaTool not initialized properly" << std::endl;
	// }

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
	      trkqual.setMVAValue(i_trkQualEntry._mvaTool->evalMVA(trkqual.values(), i_trkQualEntry._mvaMask));

	    }
	    else {
	      trkqual.setMVAStatus(MVAStatus::filled);
	    }
	  }
	  tqcol->push_back(trkqual);

	  // Get the efficiency cut that this track passes
	  Float_t passEff = -1;
	  if (_effCalib.size() > 0) {
	    for (const auto& i_pair : _effCalib) {
	      if (trkqual.MVAValue() >= i_pair.second) {
		passEff = i_pair.first;
		break;
	      }
	    }
	  }
	  rqcol->push_back(RecoQual(trkqual.status(),trkqual.MVAValue(), passEff));
	}

	if ( (tqcol->size() != rqcol->size()) || (tqcol->size() != kalSeeds.size()) ) {
	  throw cet::exception("TrackQuality") << "KalSeed, TrkQual and RecoQual sizes are inconsistent (" << kalSeeds.size() << ", " << tqcol->size() << ", " << rqcol->size() << " respectively)";
	}
      }
    }

    // put the output products into the event
    event.put(move(tqcol));
    event.put(move(rqcol));
  }
}// mu2e

DEFINE_ART_MODULE(mu2e::TrackQuality);
