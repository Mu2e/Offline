//
// Create a TrkQual object
// using TMVA::SOFIE
//
// Original author A. Edmonds
//

// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art_root_io/TFileService.h"
#include "art/Utilities/make_tool.h"
// utilities
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/Mu2eUtilities/inc/MVATools.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"
// data
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/MVAResult.hh"
#include "Offline/TrkDiag/inc/TrkQual_ANN1.hxx"
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

namespace TMVA_SOFIE_TrkQual_ANN1 {
  class Session;
}

namespace mu2e
{

  class TrackQuality : public art::EDProducer
  {
    public:
      struct Config {
        using Name=fhicl::Name;
        using Comment=fhicl::Comment;

        fhicl::Atom<art::InputTag> kalSeedPtrTag{Name("KalSeedPtrCollection"), Comment("Input tag for KalSeedPtrCollection")};
        fhicl::Atom<bool> printMVA{Name("PrintMVA"), Comment("Print the MVA used"), false};
        fhicl::Atom<std::string> datFilename{Name("datFilename"), Comment("Filename for the .dat file to use")};
        fhicl::Atom<int> debug{Name("debugLevel"), Comment("Debug printout level"), 0};
      };

      using Parameters = art::EDProducer::Table<Config>;
      TrackQuality(const Parameters& conf);

    private:
      void produce(art::Event& event) override;
      void initializeMVA(std::string xmlfilename);

      art::InputTag _kalSeedPtrTag;
      bool _printMVA;
      int _debug;

    std::shared_ptr<TMVA_SOFIE_TrkQual_ANN1::Session> mva_;

  };

  TrackQuality::TrackQuality(const Parameters& conf) :
    art::EDProducer{conf},
    _kalSeedPtrTag(conf().kalSeedPtrTag()),
    _printMVA(conf().printMVA()),
    _debug(conf().debug())
    {
      produces<MVAResultCollection>();

      ConfigFileLookupPolicy configFile;
      mva_ = std::make_shared<TMVA_SOFIE_TrkQual_ANN1::Session>(configFile(conf().datFilename()));
    }

  void TrackQuality::produce(art::Event& event ) {
    // create output
    unique_ptr<MVAResultCollection> mvacol(new MVAResultCollection());

    // get the KalSeedPtrs
    art::Handle<KalSeedPtrCollection> kalSeedPtrHandle;
    event.getByLabel(_kalSeedPtrTag, kalSeedPtrHandle);
    const auto& kalSeedPtrs = *kalSeedPtrHandle;

    // Go through the tracks and calculate their track qualities
    for (const auto& kalSeedPtr : kalSeedPtrs) {
      const auto& kalSeed = *kalSeedPtr;
      std::array<float,7> features; // the features we trained on

      // fill the hit count variables
      int nhits = 0; int nactive = 0; int ndouble = 0; int ndactive = 0; int nnullambig = 0;
      static StrawHitFlag active(StrawHitFlag::active);
        if(_debug > 1) printf("[TrackQuality::%s::%s] Printing track hit information\n", __func__, moduleDescription().moduleLabel().c_str());
      for (auto ihit = kalSeed.hits().begin(); ihit != kalSeed.hits().end(); ++ihit) {
        ++nhits;
        if (ihit->flag().hasAllProperties(active)) {
          ++nactive;
          if (ihit->ambig()==0) {
            ++nnullambig;
          }
        }
        if(_debug > 1) printf(" Hit %2i: active = %o, null = %o\n",
                              nhits, ihit->flag().hasAllProperties(active), ihit->flag().hasAllProperties(active) && ihit->ambig()==0);
        auto jhit = ihit; jhit++;
        if(jhit != kalSeed.hits().end() && ihit->strawId().uniquePanel() ==
           jhit->strawId().uniquePanel()){
          ++ndouble;
          if(ihit->flag().hasAllProperties(active)) { ++ndactive; }
        }
      }

      int ndof = nactive -5;
      if (kalSeed.hasCaloCluster()) {
        ++ndof;
      }

      int nmat = 0; int nmatactive = 0; int radlen = 0.0;
      for (std::vector<TrkStraw>::const_iterator i_straw = kalSeed.straws().begin(); i_straw != kalSeed.straws().end(); ++i_straw) {
        ++nmat;
        if (i_straw->active()) {
          ++nmatactive;
          radlen += i_straw->_radlen;
        }
      }

      features[0] = nactive;
      features[1] = (double) nactive / nhits;
      features[3] = (double) nnullambig / nactive;
      features[4] = kalSeed.fitConsistency();
      features[6] = (double)nmatactive / nactive;

      // Now get the features that are for the entrance of the trackre
      bool entrance_found = false;
      for(size_t ikinter = 0; ikinter < kalSeed.intersections().size(); ++ikinter){
        auto const& kinter = kalSeed.intersections()[ikinter];
        if (kinter.surfaceId() == SurfaceIdDetail::TT_Front) { // we only want the tracker entrance (sid=0)
          features[2] = sqrt(kinter.loopHelix().paramVar(KinKal::LoopHelix::t0_));
          features[5] = kinter.momerr();
          entrance_found = true;
          break;
        }
      }
      if (!entrance_found) {
        features[2] = -9999;
        features[5] = -9999;
      }

      std::vector<float> mvaout = mva_->infer(features.data());

      if (!entrance_found) {
        mvaout[0] = 0; // this is not a good track
      }

      if(_debug > 0) {
        printf("[TrackQuality::%s::%s] Inputs = %.0f, %.4f, %.4f, %.4f, %.4f, %.4f %.4f --> output = %.4f\n",
               __func__, moduleDescription().moduleLabel().c_str(),
               features[0], features[1], features[2], features[3], features[4], features[5], features[6], mvaout[0]);
      }

      mvacol->push_back(MVAResult(mvaout[0]));
    }

    if ( (mvacol->size() != kalSeedPtrs.size()) ) {
      throw cet::exception("TrackQuality") << "KalSeedPtr and MVAResult sizes are inconsistent (" << kalSeedPtrs.size() << ", " << mvacol->size();
    }


    // put the output products into the event
    event.put(move(mvacol));
  }
}// mu2e

DEFINE_ART_MODULE(mu2e::TrackQuality)
