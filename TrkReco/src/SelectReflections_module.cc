//
//  Select KalSeed pairs that are candidates for a reflecting particle. Input collections must have the same particle type, and have been extrapolated to the
//  tracker entrance in both directions
//  original author: D. Brown (LBNL) 2024
//
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/DelegatedParameter.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "Offline/TrkReco/inc/KalSeedSelector.hh"
// mu2e data products
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
// C++
#include <vector>
#include <tuple>
#include <limits>

namespace mu2e {
  // sorting function
  class SortKalSeeds    {
    public:
      SortKalSeeds(KalSeedSelector const& sel): sel_(sel) {}
      bool operator()(KalSeedPtr const& a, KalSeedPtr const& b) const { return sel_.isBetter(*b,*a); }
    private:
      KalSeedSelector const& sel_;
  };


  class SelectReflections : public art::EDProducer {
    public:
      enum bestpair {mom=0, deltat,deltap}; // selection options for defining the best pair
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
        fhicl::Atom<int> debug{ Name("debugLevel"), Comment("Debug Level"), 0};
        fhicl::Atom<art::InputTag> upstreamTag {Name("UpstreamKalSeedCollection"), Comment("Upstream KalSeed collection") };
        fhicl::Atom<art::InputTag> downstreamTag {Name("DownstreamKalSeedCollection"), Comment("Downstream KalSeed collection") };
        fhicl::Atom<double> maxdt{ Name("MaxDeltaT"), Comment("Maximum time difference at tracker entrance")};
        fhicl::Atom<double> maxdp{ Name("MaxDeltaP"), Comment("Maximum scalar momentum difference at tracker entrance")};
        fhicl::DelegatedParameter selector{ Name("Selector"), Comment("Selector parameters")};
        fhicl::Atom<int> selectbest{ Name("SelectBestPair"), Comment("Method to select the best pair if multiple pairs are found ")};
      };
      using Parameters = art::EDProducer::Table<Config>;
      explicit SelectReflections(const Parameters& config);
      void produce(art::Event& evt) override;
    private:
      int debug_;
      double maxdt_, maxdp_;
      art::ProductToken<KalSeedCollection> upstreamtoken_, downstreamtoken_;
      std::unique_ptr<KalSeedSelector> selector_;
      bestpair selbest_;
  };

  SelectReflections::SelectReflections(const Parameters& config) : art::EDProducer{config},
    debug_(config().debug()),
    maxdt_(config().maxdt()),
    maxdp_(config().maxdp()),
    upstreamtoken_ { consumes<KalSeedCollection> (config().upstreamTag()) },
    downstreamtoken_ { consumes<KalSeedCollection> (config().downstreamTag()) },
    selector_(art::make_tool<KalSeedSelector>(config().selector.get<fhicl::ParameterSet>())),
    selbest_(static_cast<bestpair>(config().selectbest())) {
      produces<KalSeedPtrCollection> ();
    }

  void SelectReflections::produce(art::Event& event) {
    static const SurfaceId trkfront (SurfaceIdDetail::TT_Front);
    // create output
    std::unique_ptr<KalSeedPtrCollection> mkseeds(new KalSeedPtrCollection);
    auto const& upksch = event.getValidHandle<KalSeedCollection>(upstreamtoken_);
    auto const& upksc = *upksch;
    auto const& downksch = event.getValidHandle<KalSeedCollection>(downstreamtoken_);
    auto const& downksc = *downksch;
    if(debug_ > 0) std::cout << "Upstream " << upksc.size() << " , Downstream " << downksc.size() << " KalSeeds" << std::endl;
    if(upksc.size() > 0 && downksc.size() > 0){
      std::vector<std::tuple<size_t, size_t,double, double, double>> matches; // matched up and downstream track, with (downstream) mom, dt and dmom
      for(size_t iup = 0; iup <upksc.size(); ++iup){
        auto const& upks = upksc[iup];
        if(selector_->select(upks)){
          // find the appropriate intersection for comparison
          auto uptrkiinter = upks.intersections().end();
          for(auto upiinter = upks.intersections().begin(); upiinter != upks.intersections().end(); ++upiinter){
            auto const& upinter = *upiinter;
            if(upinter.surfaceId() == trkfront && upinter.momentum3().Z() > 0.0){ // correct surface and direction
              if(debug_ > 1) std::cout << "Found upstream intersection mom " << upinter.mom() << " time " << upinter.time() << std::endl;
              uptrkiinter = upiinter;
              break;
            }
          }
          // if no intersections found, skip testing for a match with this track
          if(uptrkiinter == upks.intersections().end())break;
          // otherwise, search for a matching downstream track
          for(size_t idown = 0; idown <downksc.size(); ++idown){
            auto const& downks = downksc[idown];
            if(selector_->select(downks)){
              // find the appropriate intersection for comparison
              auto downtrkiinter = downks.intersections().end();
              for(auto downiinter = downks.intersections().begin(); downiinter != downks.intersections().end(); ++downiinter){
                auto const& downinter = *downiinter;
                if(downinter.surfaceId() == trkfront && downinter.momentum3().Z() < 0.0){ // correct surface and direction
                  if(debug_ > 1) std::cout << "Found downstream intersection mom " << downinter.mom() << " time " << downinter.time() << std::endl;
                  downtrkiinter = downiinter;
                  break;
                }
              }
              if(downtrkiinter != downks.intersections().end()){
                // potentially matching tracks: compare time and momentum
                double dt = fabs(uptrkiinter->time() - downtrkiinter->time());
                double dmom = fabs(uptrkiinter->mom() - downtrkiinter->mom());
                if( dt < maxdt_ && dmom < maxdp_){
                  if(debug_ > 1) std::cout << "Found matching track pair " << std::endl;
                  matches.emplace_back(iup,idown,downtrkiinter->mom(),dt,dmom);
                }
              }
            }
          }
        }
      }
      // if there are >1 matches select the best
      int ibest = -1;
      if(matches.size() >0 ){
        ibest = 0;
        if(matches.size()>1){
          if(debug_ > 1) std::cout << "Selecting best reflection pair from " << matches.size() << " candidates " << std::endl;
          double value = std::numeric_limits<float>::max();
          for (size_t imatch = 0; imatch < matches.size(); ++imatch) {
            auto const& match = matches[imatch];
            if(selbest_ == mom && -std::get<2>(match) < value){
              ibest = imatch;
              value = -std::get<2>(match); // sign to make highest momentum best
            } else if(selbest_ == deltat && std::get<3>(match) < value){
              ibest = imatch;
              value = std::get<3>(match);
            } else if(selbest_ == deltap && std::get<4>(match) < value){
              ibest = imatch;
              value = std::get<4>(match);
            }
          }
        }
      }
      if(ibest > 0){
        if(debug_ > 0) std::cout << "Found Reflecting particle candidate, downstream momentum " << std::get<2>(matches[ibest])
          << " delta t " << std::get<3>(matches[ibest])
            << " delta P " << std::get<4>(matches[ibest]) << std::endl;
        mkseeds->emplace_back(upksch,std::get<0>(matches[ibest])); // store the upstream track first by convention
        mkseeds->emplace_back(downksch,std::get<1>(matches[ibest]));
      }

    }
    event.put(std::move(mkseeds));
  }
}
DEFINE_ART_MODULE(mu2e::SelectReflections)
