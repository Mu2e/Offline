//
//  Select KalSeed pairs that are candidates for being 'the same' track
//  original author: D. Brown (LBNL) 2025
//
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/DelegatedParameter.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDFilter.h"
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


  class SelectSameTrack : public art::EDFilter {
    public:
      enum bestpair {mom=0, deltat,deltap,hitfrac}; // selection options for defining the best pair
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
        fhicl::Atom<int> debug{ Name("debugLevel"), Comment("Debug Level"), 0};
        fhicl::Atom<art::InputTag> primaryTag {Name("PrimaryKalSeedCollection"), Comment("Primary KalSeed collection") };
        fhicl::Atom<art::InputTag> secondaryTag {Name("SecondaryKalSeedCollection"), Comment("Secondary KalSeed collection") };
        fhicl::Atom<double> maxdt{ Name("MaxDeltaT"), Comment("Maximum time difference at tracker mid")};
        fhicl::Atom<double> maxdp{ Name("MaxDeltaP"), Comment("Maximum scalar momentum difference at tracker mid")};
        fhicl::Atom<double> minhitfrac{ Name("MinHitFrac"), Comment("Minimum number of common tracker hits")};
        fhicl::DelegatedParameter selector{ Name("Selector"), Comment("Selector parameters")};
        fhicl::Atom<int> selectbest{ Name("SelectBestPair"), Comment("Method to select the best pair if multiple pairs are found ")};
        fhicl::Atom<std::string> comparisonDirection{ Name("ComparisonDirection"), Comment("Direction (dZ/dt) to compare fits, 'Unknown' means test all")};
        fhicl::Atom<std::string> surface{ Name("ComparisonSurface"), Comment("Surface at which to compare fits")};
      };

      struct TrkMatch {
        size_t priIndex_; // indices of matched tracks
        size_t secIndex_;
        double primom_; // primary momentum
        double dt_; // primary - secondary time difference
        double dmom_; // primary - secondary momentum difference
        double hfrac_; // matched hit fraction
        TrkMatch(size_t ipri, size_t isec, double primom, double dt, double dmom, double hfrac) :
          priIndex_(ipri), secIndex_(isec), primom_(primom), dt_(dt), dmom_(dmom), hfrac_(hfrac) {}
      };

      using Parameters = art::EDFilter::Table<Config>;
      explicit SelectSameTrack(const Parameters& config);
      bool filter(art::Event& evt) override;
    private:
      int debug_;
      double maxdt_, maxdp_, minhf_;
      SurfaceId compsurf_;
      TrkFitDirection compdir_;
      art::ProductToken<KalSeedCollection> primarytoken_, secondarytoken_;
      std::unique_ptr<KalSeedSelector> selector_;
      bestpair selbest_;
  };

  SelectSameTrack::SelectSameTrack(const Parameters& config) : art::EDFilter{config},
    debug_(config().debug()),
    maxdt_(config().maxdt()),
    maxdp_(config().maxdp()),
    minhf_(config().minhitfrac()),
    compsurf_(config().surface()),
    compdir_(config().comparisonDirection()),
    primarytoken_ { consumes<KalSeedCollection> (config().primaryTag()) },
    secondarytoken_ { consumes<KalSeedCollection> (config().secondaryTag()) },
    selector_(art::make_tool<KalSeedSelector>(config().selector.get<fhicl::ParameterSet>())),
    selbest_(static_cast<bestpair>(config().selectbest())) {
      produces<KalSeedPtrCollection> ();
    }

  bool SelectSameTrack::filter(art::Event& event) {
    // create output
    std::unique_ptr<KalSeedPtrCollection> mkseeds(new KalSeedPtrCollection);
    auto const& priksch = event.getValidHandle<KalSeedCollection>(primarytoken_);
    auto const& priksc = *priksch;
    auto const& secksch = event.getValidHandle<KalSeedCollection>(secondarytoken_);
    auto const& secksc = *secksch;
    if(debug_ > 1) std::cout << "Primary " << priksc.size() << " , Secondary " << secksc.size() << " KalSeeds" << std::endl;
    bool keep(false);
    if(priksc.size() > 0 && secksc.size() > 0){
      std::vector<TrkMatch> matches;
      for(size_t ipri = 0; ipri <priksc.size(); ++ipri){
        auto const& priks = priksc[ipri];
        if(selector_->select(priks)){
          if(debug_ > 2)std::cout << "Selected primary track " << std::endl;
          // find the appropriate intersection for comparison
          auto pritrkiinter = priks.intersections().end();
          for(auto priiinter = priks.intersections().begin(); priiinter != priks.intersections().end(); ++priiinter){
            auto const& priinter = *priiinter;
            if(priinter.surfaceId() == compsurf_ && priinter.momentum3().Z()*compdir_.dzdt() > 0.0){ // correct surface and direction
              if(debug_ > 1) std::cout << "Found primary intersection mom " << priinter.momentum3() << " time " << priinter.time() << std::endl;
              pritrkiinter = priiinter;
              break;
            }
          }
          // if no intersections found, skip testing for a match with this track
          if(pritrkiinter == priks.intersections().end())break;
          // otherwise, search for a matching secondary track
          for(size_t isec = 0; isec <secksc.size(); ++isec){
            auto const& secks = secksc[isec];
            if(selector_->select(secks)){
              if(debug_ > 2)std::cout << "Selected secondary track " << std::endl;
              // find the appropriate intersection for comparison
              auto sectrkiinter = secks.intersections().end();
              for(auto seciinter = secks.intersections().begin(); seciinter != secks.intersections().end(); ++seciinter){
                auto const& secinter = *seciinter;
                if(secinter.surfaceId() == compsurf_ && secinter.momentum3().Z()*compdir_.dzdt() > 0.0){ // correct surface and direction
                  if(debug_ > 1) std::cout << "Found secondary intersection mom " << secinter.momentum3() << " time " << secinter.time() << std::endl;
                  sectrkiinter = seciinter;
                  break;
                }
              }
              if(sectrkiinter != secks.intersections().end()){
                // potentially matching tracks: compare time and momentum
                double dt = fabs(pritrkiinter->time() - sectrkiinter->time());
                double dmom = fabs(pritrkiinter->mom() - sectrkiinter->mom());
                if( dt < maxdt_ && dmom < maxdp_){
                  // now test hit overlap
                  std::set<unsigned> prihits;
                  for(auto const& hit: priks.hits()){
                    if(hit._ambig > WireHitState::inactive)prihits.insert(hit._index);
                  }
                  unsigned nover(0);
                  for(auto const& hit: secks.hits()){
                    if(hit._ambig > WireHitState::inactive && prihits.find(hit._index) != prihits.end())++nover;
                  }
                  double hfrac = float(nover)/float(prihits.size());
                  if(hfrac > minhf_){
                    if(debug_ > 1) std::cout << "Found matching track pair " << std::endl;
                    matches.emplace_back(ipri,isec,pritrkiinter->mom(),dt,dmom,hfrac);
                  }
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
            if(selbest_ == mom && -match.primom_ < value){
              ibest = imatch;
              value = -match.primom_; // sign to make highest momentum best
            } else if(selbest_ == deltat && match.dt_ < value){
              ibest = imatch;
              value = match.dt_;
            } else if(selbest_ == deltap && match.dmom_ < value){
              ibest = imatch;
              value = match.dmom_;
            } else if(selbest_ == hitfrac && match.hfrac_ < value){
              ibest = imatch;
              value = match.hfrac_;
            }
          }
        }
      }
      if(ibest > -1){
        if(debug_ > 0) std::cout << "Found track pair candidate, primary momentum " << matches[ibest].primom_
          << " delta t " << matches[ibest].dt_
          << " delta P " << matches[ibest].dmom_
          << " hitfrac " << matches[ibest].hfrac_
          << std::endl;
        mkseeds->emplace_back(priksch,matches[ibest].priIndex_); // store the primary track first by convention
        mkseeds->emplace_back(secksch,matches[ibest].secIndex_);
        keep = true;
      }
    }
    event.put(std::move(mkseeds));
    return keep;
  }
}
DEFINE_ART_MODULE(mu2e::SelectSameTrack)
