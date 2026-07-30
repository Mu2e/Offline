//
//  Select KalSeed pairs that are candidates for a reflecting particle. Input collections must have the same particle type, and have been extrapolated to the
//  tracker entrance in both directions
//  original author: D. Brown (LBNL) 2024
//
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
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


  class SelectReflections : public art::EDFilter {
    public:
      enum bestpair {mom=0, deltat,deltap, nactive}; // selection options for defining the best pair
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
        fhicl::Atom<int> debug{ Name("debugLevel"), Comment("Debug Level"), 0};
        fhicl::Atom<art::InputTag> upstreamTag {Name("UpstreamKalSeedCollection"), Comment("Upstream KalSeed collection") };
        fhicl::Atom<art::InputTag> dnstreamTag {Name("DownstreamKalSeedCollection"), Comment("Downstream KalSeed collection") };
        fhicl::OptionalAtom<std::string> surface{ Name("Surface"), Comment("Surface to compare fits at")};
        fhicl::OptionalAtom<bool> ptrs{ Name("PtrCollection"), Comment("Inputs are KalSeedPtr collections")};
        fhicl::Atom<double> maxdt{ Name("MaxDeltaT"), Comment("Maximum time difference at comparison surface")};
        fhicl::Atom<double> maxdp{ Name("MaxDeltaP"), Comment("Maximum scalar momentum difference at comparison surface")};
        fhicl::DelegatedParameter selector{ Name("Selector"), Comment("Selector parameters")};
        fhicl::Atom<int> selectbest{ Name("SelectBestPair"), Comment("Method to select the best pair if multiple pairs are found ")};
      };
      using Parameters = art::EDFilter::Table<Config>;
      explicit SelectReflections(const Parameters& config);
      bool filter(art::Event& evt) override;
      void endJob() override;
    private:
      int debug_;
      art::InputTag upstreamTag_, dnstreamTag_;
      SurfaceId sid_;
      bool sidmatch_;
      bool ptrs_;
      double maxdt_, maxdp_;
      std::unique_ptr<KalSeedSelector> selector_;
      bestpair selbest_;
      unsigned nevts_=0, nref_=0, nmultref_=0;
  };

  SelectReflections::SelectReflections(const Parameters& config) : art::EDFilter{config},
    debug_(config().debug()),
    upstreamTag_(config().upstreamTag()),
    dnstreamTag_(config().dnstreamTag()),
    ptrs_(config().ptrs()),
    maxdt_(config().maxdt()),
    maxdp_(config().maxdp()),
    selector_(art::make_tool<KalSeedSelector>(config().selector.get<fhicl::ParameterSet>())),
    selbest_(static_cast<bestpair>(config().selectbest()))
    {
      produces<KalSeedPtrCollection> ();
      if(ptrs_){
        consumes<KalSeedPtrCollection> (upstreamTag_);
        consumes<KalSeedPtrCollection> (dnstreamTag_);
      } else {
        consumes<KalSeedCollection> (upstreamTag_);
        consumes<KalSeedCollection> (dnstreamTag_);
      }
      std::string sidname;
      sidmatch_ = config().surface(sidname);
      if(sidmatch_)sid_ = SurfaceId(sidname);
    }

  bool SelectReflections::filter(art::Event& event) {
    ++nevts_;
    const static SurfaceId usid;
    // create output
    std::unique_ptr<KalSeedPtrCollection> mkseeds(new KalSeedPtrCollection);

    // input
    KalSeedPtrCollection upksptrs, dnksptrs;
    art::Handle<KalSeedCollection> upksh,dnksh;;

    if(ptrs_){
      auto upksptrsh = event.getValidHandle<KalSeedPtrCollection>(upstreamTag_);
      upksptrs = *upksptrsh;
      auto dnksptrsh = event.getValidHandle<KalSeedPtrCollection>(dnstreamTag_);
      dnksptrs = *dnksptrsh;
    } else {
      upksh = event.getHandle<KalSeedCollection>(upstreamTag_);
      for(size_t iks = 0; iks < upksh->size(); ++iks){
        upksptrs.emplace_back(upksh,iks);
      }
      dnksh = event.getHandle<KalSeedCollection>(dnstreamTag_);
      for(size_t iks = 0; iks < dnksh->size(); ++iks){
        dnksptrs.emplace_back(dnksh,iks);
      }
    }

    if(debug_ > 1) std::cout << "Upstream " << upksptrs.size() << " , Downstream " << dnksptrs.size() << " KalSeeds" << std::endl;
    bool keep(false);
    if(upksptrs.size() > 0 && dnksptrs.size() > 0){
      std::vector<std::tuple<size_t, size_t,double, double, double, unsigned>> matches; // matched up and dnstream track, with (dnstream) mom, dt and dmom
      for(size_t iup = 0; iup <upksptrs.size(); ++iup){
        auto const& upks = *upksptrs[iup];
        double upt0;
        auto upmom = upks.t0Segment(upt0)->momentum3();
        if(upmom.Z() < 0 && selector_->select(upks)){
          if(debug_ > 2)std::cout << "Selected upstream track " << std::endl;
          auto upinters = upks.intersections(sid_);
          // search for a matching dnstream track
          for(size_t idn = 0; idn <dnksptrs.size(); ++idn){
            auto const& dnks = *dnksptrs[idn];
            double dnt0;
            auto dnmom = dnks.t0Segment(dnt0)->momentum3();
            if(dnmom.Z() > 0 && selector_->select(dnks)){
              if(debug_ > 2)std::cout << "Selected dnstream track " << std::endl;
              if(sidmatch_) {
                // find the appropriate intersection for comparison
                auto dninters = dnks.intersections(sid_);
                for(auto upinter : upinters){
                  auto upimom = upinter->momentum3();
                  for(auto dninter : dninters){
                    auto dnimom = dninter->momentum3();
                    if(upimom.Z()*dnimom.Z() > 0.0){ // same direction
                      if(debug_ > 1) std::cout << "Same-direction upstream intersection mom " << upimom << " time " << upinter->time()
                        << " dnstream intersection mom " << dnimom << " time " << dninter->time() << std::endl;
                      // potentially matching tracks: compare time and momentum
                      double dt = fabs(upinter->time() - dninter->time());
                      double dmom = fabs(upimom.R() - dnimom.R());
                      if(debug_ > 2) std::cout << "Testing intersection track pair, dt " << dt << " dmom " << dmom << std::endl;
                      if( dt < maxdt_ && dmom < maxdp_){
                        if(debug_ > 1) std::cout << "Found matching track pair, dt " << dt << " dmom " << dmom << std::endl;
                        matches.emplace_back(iup,idn,dnimom.R(),dt,dmom,upks.nHits()+dnks.nHits());
                      }
                    }
                  }
                }
              } else {
                // match just using t0
                double dt = fabs(upt0 - dnt0);
                double dmom = fabs(dnmom.R() - upmom.R());
                if(debug_ > 2) std::cout << "Testing t0 track pair, dt " << dt << " dmom " << dmom << std::endl;
                if( dt < maxdt_ && dmom < maxdp_){
                  if(debug_ > 1) std::cout << "Found matching track pair, dt " << dt << " dmom " << dmom << std::endl;
                  matches.emplace_back(iup,idn,dnmom.R(),dt,dmom,upks.nHits()+dnks.nHits());
                }
              }
            }
          }
        }
      }
      // if there are >1 matches select the best
      int ibest = -1;
      if(matches.size() >0 ){
        ++nref_;
        ibest = 0;
        if(matches.size()>1){
          ++nmultref_;
          if(debug_ > 0) std::cout << "Selecting best reflection pair from " << matches.size() << " candidates " << " according to " << selbest_ << " event " << event.id() << std::endl;
          double value = (selbest_ == deltat || selbest_ == deltap) ? std::numeric_limits<double>::max() : 0;
          for (size_t imatch = 0; imatch < matches.size(); ++imatch) {
            auto const& match = matches[imatch];
            if(debug_ > 1)std::cout << "Match " << imatch << " has dnstream momentum " << std::get<2>(match) << " dt " << std::get<3>(match) << " dp " << std::get<4>(match) << " nactive " << std::get<5>(match) << std::endl;
            if(selbest_ == mom && std::get<2>(match) > value){
              ibest = imatch;
              value = std::get<2>(match);
            } else if(selbest_ == deltat && std::get<3>(match) < value){
              ibest = imatch;
              value = std::get<3>(match);
            } else if(selbest_ == deltap && std::get<4>(match) < value){
              ibest = imatch;
              value = std::get<4>(match);
            } else if(selbest_ == nactive && std::get<5>(match) > value){
              ibest = imatch;
              value = std::get<5>(match);
            }
          }
        }
      }
      if(ibest > -1){
        if(debug_ > 0) std::cout << "Found Reflecting particle candidate, dnstream momentum " << std::get<2>(matches[ibest])
          << " delta t " << std::get<3>(matches[ibest])
            << " delta P " << std::get<4>(matches[ibest]) << std::endl;

        if(ptrs_){
          mkseeds->emplace_back(upksptrs[std::get<0>(matches[ibest])]); // store the upstream track first by convention
          mkseeds->emplace_back(dnksptrs[std::get<1>(matches[ibest])]);
        } else {
          mkseeds->emplace_back(upksh,std::get<0>(matches[ibest])); // store the upstream track first by convention
          mkseeds->emplace_back(dnksh,std::get<1>(matches[ibest]));
        }
        keep = true;
      }
    }
    event.put(std::move(mkseeds));
    return keep;
  }

  void SelectReflections::endJob() {
    std::cout << moduleDescription().moduleLabel() << " processed " << nevts_ << " events and found " << nref_ << " reflections with " << nmultref_ << " multiple reflections" << std::endl;
  }
}
DEFINE_ART_MODULE(mu2e::SelectReflections)
