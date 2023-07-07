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
#include "art_root_io/TFileService.h"
#include "art/Utilities/make_tool.h"
// utilities
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/Mu2eUtilities/inc/MVATools.hh"
#include "Offline/AnalysisConditions/inc/TrkQualCatalog.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
// data
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/TrkQual.hh"
#include "Offline/RecoDataProducts/inc/RecoQual.hh"
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
    auto const& ptable = GlobalConstantsHandle<ParticleDataList>();
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
      trkQualEntry._mvaTool.showMVA();
    }

    // Go through the tracks and calculate their track qualities
    for (const auto& i_kalSeed : kalSeeds) {
      TrkQual trkqual;

      static TrkFitFlag goodfit(TrkFitFlag::kalmanOK);
      if (i_kalSeed.status().hasAllProperties(goodfit)) {

        // fill the hit count variables
        int nhits = 0; int nactive = 0; int ndouble = 0; int ndactive = 0; int nnullambig = 0;
        static StrawHitFlag active(StrawHitFlag::active);
        for (auto ihit = i_kalSeed.hits().begin(); ihit != i_kalSeed.hits().end(); ++ihit) {
          ++nhits;
          if (ihit->flag().hasAllProperties(active)) {
            ++nactive;
            if (ihit->ambig()==0) {
              ++nnullambig;
            }
          }
          auto jhit = ihit; jhit++;
          if(jhit != i_kalSeed.hits().end() && ihit->strawId().uniquePanel() ==
              jhit->strawId().uniquePanel()){
            ++ndouble;
            if(ihit->flag().hasAllProperties(active)) { ++ndactive; }
          }
        }

        int ndof = nactive -5;
        if (i_kalSeed.hasCaloCluster()) {
          ++ndof;
        }

        int nmat = 0; int nmatactive = 0; int radlen = 0.0;
        for (std::vector<TrkStraw>::const_iterator i_straw = i_kalSeed.straws().begin(); i_straw != i_kalSeed.straws().end(); ++i_straw) {
          ++nmat;
          if (i_straw->active()) {
            ++nmatactive;
            radlen += i_straw->radLen();
          }
        }


        trkqual[TrkQual::nactive] = nactive;
        trkqual[TrkQual::factive] = (double) nactive / nhits;
        trkqual[TrkQual::fdouble] = (double) ndactive / nactive;
        trkqual[TrkQual::fnullambig] = (double) nnullambig / nactive;
        trkqual[TrkQual::fstraws] = (double)nmatactive / nactive;

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
          double charge = ptable->particle(i_kalSeed.particle()).charge();
          trkqual[TrkQual::momerr] = bestkseg->momerr();
          trkqual[TrkQual::d0] = -1*charge*bestkseg->helix().d0();
          trkqual[TrkQual::rmax] = -1*charge*(bestkseg->helix().d0() + 2.0/bestkseg->helix().omega());

          trkqual.setMVAStatus(MVAStatus::calculated);
          trkqual.setMVAValue(trkQualEntry._mvaTool.evalMVA(trkqual.values(), trkQualEntry._mvaMask));

        }
        else {
          trkqual.setMVAStatus(MVAStatus::filled);
        }
      }
      tqcol->push_back(trkqual);

      // Get the efficiency cut that this track passes
      Float_t passCalib = 0.0; // everything will pass a 100% efficient cut
      if (trkQualEntry._calibrated) {
        //      std::cout << _trainName << " = " << trkqual.MVAValue() << std::endl;
        for (const auto& i_pair : trkQualEntry._effCalib) {
          //      std::cout << i_pair.first << ", " << i_pair.second << ": ";
          if (trkqual.MVAValue() >= i_pair.second) {
            //      std::cout << "PASSES" << std::endl;
            passCalib = i_pair.first;
          }
          else {
            //      std::cout << "FAILS" << std::endl;
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

DEFINE_ART_MODULE(mu2e::TrackQuality)
