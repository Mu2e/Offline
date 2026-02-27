// "Cluster" calo detector steps together and then filter on these clusters
// Steps are first clustered in time using a moving window
// Time clustered hits are then histogrammed in (x,y) to search for spacial clusters
// Original author; Michael MacKenzie, 2026

// Framework
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/OptionalTable.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// cet
#include "cetlib_except/exception.h"

// Offline
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/CalorimeterGeom/inc/CaloGeomUtil.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"
#include "Offline/MCDataProducts/inc/CaloShowerStep.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/PrimaryParticle.hh"

// C++
#include <string>
#include <map>
#include <sstream>

// ROOT
#include "TH2.h"

namespace mu2e {

  class CaloDtsClusterFilter : public art::EDFilter {
  public:

    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::Sequence<art::InputTag> caloSteps         { Name("CaloShowerSteps")      , Comment("CaloShowerStep collections") };
      fhicl::Atom<double>            minCaloStepEnergy { Name("MinimumStepEnergy")    , Comment("Minimum Calo step energy"), 0.0 };
      fhicl::Atom<double>            minCaloStepTime   { Name("MinimumStepTime")      , Comment("Minimum Calo step time")};
      fhicl::Atom<double>            maxCaloStepTime   { Name("MaximumStepTime")      , Comment("Maximum Calo step time")};
      fhicl::Atom<double>            minClusterEnergy  { Name("MinimumClusterEnergy") , Comment("Minimum cluster energy")};
      fhicl::Atom<double>            timeWindow        { Name("TimeWindow")           , Comment("Time window for clustering")};
      fhicl::Atom<double>            spaceWindow       { Name("SpaceWindow")          , Comment("Space window for clustering")};
      fhicl::Atom<bool>              nullFilter        { Name("NullFilter")           , Comment("If true, do not apply any filtering and accept all events"), false };
      fhicl::Atom<int>               diagLevel         { Name("DiagLevel")            , Comment("Diagnostic output level"), 0 };
    };

    using Parameters = art::EDFilter::Table<Config>;
    explicit CaloDtsClusterFilter(const Parameters& conf);
    bool filter(art::Event& event) override;
    bool beginRun(art::Run& aRun) override;
    void endJob() override;

  private:
    bool goodParticle(SimParticle const& simp) const;
    CLHEP::Hep3Vector stepPosition(const CaloShowerStep& step) const;
    void sortStepsByTime(std::vector<const CaloShowerStep*>& steps) const;
    void fillStepsList(std::vector<const CaloShowerStep*>& steps, const art::Event& event) const;
    bool checkClusterHistogram(const std::unique_ptr<TH2D>& hist) const;
    bool findSpacialCluster(const std::vector<const CaloShowerStep*>& steps, size_t first_index, size_t last_index);
    bool findCluster(const std::vector<const CaloShowerStep*>& steps);

    // Input data
    std::vector<art::InputTag> caloStepsTags_;
    double                     minCaloStepEnergy_;
    double                     minCaloStepTime_;
    double                     maxCaloStepTime_;
    double                     minClusterEnergy_;
    double                     timeWindow_;
    double                     spaceWindow_;
    bool                       nullFilter_;
    int                        diagLevel_;

    // Data
    const Calorimeter*         calorimeter_;
    std::unique_ptr<TH2D>      clusterHistD0_;
    std::unique_ptr<TH2D>      clusterHistD1_;
    long unsigned              nEvents_;
    long unsigned              nPassed_;
  };

  //================================================================
  CaloDtsClusterFilter::CaloDtsClusterFilter(const Parameters& conf)
    : art::EDFilter{conf}
    , caloStepsTags_(conf().caloSteps())
    , minCaloStepEnergy_(conf().minCaloStepEnergy())
    , minCaloStepTime_(conf().minCaloStepTime())
    , maxCaloStepTime_(conf().maxCaloStepTime())
    , minClusterEnergy_(conf().minClusterEnergy())
    , timeWindow_(conf().timeWindow())
    , spaceWindow_(conf().spaceWindow())
    , nullFilter_(conf().nullFilter())
    , diagLevel_(conf().diagLevel())
    , calorimeter_(nullptr)
    , nEvents_(0)
    , nPassed_(0)
  {
    if(caloStepsTags_.size() == 0)          throw cet::exception("BADCONFIG") << "At least one CaloShowerStep collection must be specified\n";
    if(minCaloStepTime_ < 0.)               throw cet::exception("BADCONFIG") << "Minimum Calo step time must be >= 0\n";
    if(maxCaloStepTime_ < 0.)               throw cet::exception("BADCONFIG") << "Maximum Calo step time must be >= 0\n";
    if(minCaloStepTime_ > maxCaloStepTime_) throw cet::exception("BADCONFIG") << "Minimum Calo step time must be <= Maximum Calo step time\n";
    if(minClusterEnergy_ < 0.)              throw cet::exception("BADCONFIG") << "Minimum cluster energy must be >= 0\n";
    if(timeWindow_ < 0.)                    throw cet::exception("BADCONFIG") << "Time window must be >= 0\n";
    if(spaceWindow_ < 0.)                   throw cet::exception("BADCONFIG") << "Space window must be >= 0\n";

    constexpr double hist_range = 1000.; // -1000 - 1000 mm in (x,y)
    const int n_space_bins = 2.*hist_range/spaceWindow_; // bin width = space window
    clusterHistD0_ = std::make_unique<TH2D>("clusterHistD0", "Hit x vs. y",
                                          n_space_bins, -hist_range, hist_range,
                                          n_space_bins, -hist_range, hist_range);
    clusterHistD1_ = std::make_unique<TH2D>("clusterHistD1", "Hit x vs. y",
                                          n_space_bins, -hist_range, hist_range,
                                          n_space_bins, -hist_range, hist_range);
  }

  //--------------------------------------------------------------------------------
  bool CaloDtsClusterFilter::goodParticle(SimParticle const& simp) const {
    // Check if the SimParticle is a valid particle to consider for pileup filtering
    return true; // For now, accept all particles
  }

  //--------------------------------------------------------------------------------
  CLHEP::Hep3Vector CaloDtsClusterFilter::stepPosition(const CaloShowerStep& step) const {
    // Get the position of the CaloShowerStep in tracker coordinates, so (x,y) are centered on (0,0)
    const CLHEP::Hep3Vector& pos = step.position();
    const auto& calo_geom = calorimeter_->geomUtil();
    const auto pos_in_mu2e = calo_geom.crystalToMu2e(step.volumeG4ID(), pos);
    const auto pos_in_trk = calo_geom.mu2eToTracker(pos_in_mu2e);
    return pos_in_trk;
  }

  //--------------------------------------------------------------------------------
  void CaloDtsClusterFilter::sortStepsByTime(std::vector<const CaloShowerStep*>& steps) const {
    const double mbtime = GlobalConstantsHandle<PhysicsParams>()->getNominalDRPeriod();
    std::sort(steps.begin(), steps.end(), [mbtime](const CaloShowerStep* a, const CaloShowerStep* b) {
      return std::fmod(a->time(), mbtime) < std::fmod(b->time(), mbtime);
    });
  }

  //--------------------------------------------------------------------------------
  void CaloDtsClusterFilter::fillStepsList(std::vector<const CaloShowerStep*>& steps, const art::Event& event) const {
    steps.clear();
    const double mbtime = GlobalConstantsHandle<PhysicsParams>()->getNominalDRPeriod();
    for(const auto& tag : caloStepsTags_) {
      const auto handle = event.getValidHandle<CaloShowerStepCollection>(tag);
      for(const auto& css : *handle) {
        const auto sim_step = css.simParticle();
        if(sim_step.isNull()) continue;
        if(!goodParticle(*sim_step)) continue; // skip steps from particles we don't care about
        if(css.energyDepBirks() < minCaloStepEnergy_) continue; // skip low energy steps
        double wrapped_time = std::fmod(css.time(), mbtime);
        if(wrapped_time < 0.) wrapped_time += mbtime; // ensure time is in [0, mbtime)
        if(wrapped_time < minCaloStepTime_) continue; // skip steps that are too early in time
        if(css.time() > maxCaloStepTime_) continue; // skip steps that are too far in time to properly apply fmod
        steps.push_back(&css);
      }
    }
  }

  //--------------------------------------------------------------------------------
  bool CaloDtsClusterFilter::checkClusterHistogram(const std::unique_ptr<TH2D>& hist) const {
    // Quickly check if there are any bins above the minimum
    const double hist_max = hist->GetMaximum();
    if(hist_max > minClusterEnergy_) return true; // found a bin above the threshold, accept
    if(hist_max < minClusterEnergy_/2.) return false; // no bins even close to the threshold, reject

    // Loop over the histogram bins and check if there are any bins above the energy threshold
    const int n_bins_x = hist->GetNbinsX();
    const int n_bins_y = hist->GetNbinsY();
    for(int x_bin = 1; x_bin <= n_bins_x; ++x_bin) {
      for(int y_bin = 1; y_bin <= n_bins_y; ++y_bin) {
        const double e_bin = hist->GetBinContent(x_bin, y_bin);
        // if(e_bin < minClusterEnergy_/2.) continue; // skip bins that are too low to be interesting
        // Check if there are neighboring bins that together with this bin would pass the threshold
        for(int dx = -1; dx <= 0; ++dx) { // only need to check half the neighbors to avoid double counting
          for(int dy = -1; dy <= 1; ++dy) {
            if(dx == 0 && dy >= 0) continue; // skip the central and upper middle bin
            const int x_bin_n = x_bin + dx;
            const int y_bin_n = y_bin + dy;
            if(x_bin_n < 1 || x_bin_n > n_bins_x || y_bin_n < 1 || y_bin_n > n_bins_y) continue; // skip out of range bins
            const double e_bin_n = hist->GetBinContent(x_bin_n, y_bin_n);
            if(e_bin + e_bin_n > minClusterEnergy_) {
              if(diagLevel_ > 0) {
                std::cout << "  Found neighboring bins above threshold: (" << hist->GetXaxis()->GetBinCenter(x_bin) << ", " << hist->GetYaxis()->GetBinCenter(y_bin) << ") with e = " << e_bin
                          << " and (" << hist->GetXaxis()->GetBinCenter(x_bin_n) << ", " << hist->GetYaxis()->GetBinCenter(y_bin_n) << ") with e = " << e_bin_n
                          << " total e = " << e_bin + e_bin_n
                          << std::endl;
              }
              return true;
            }
          }
        }
      }
    }
    return false;
  }

  //--------------------------------------------------------------------------------
  bool CaloDtsClusterFilter::findSpacialCluster(const std::vector<const CaloShowerStep*>& steps,
                                                size_t first_index, size_t last_index) {

    // Fill the histograms with the step positions + energy
    clusterHistD0_->Reset();
    clusterHistD1_->Reset();
    for(size_t index = first_index; index <= last_index; ++index) {
      const auto* step_i = steps[index];
      const auto pos = stepPosition(*step_i);
      const bool disk_0 = step_i->volumeG4ID() < int(calorimeter_->nCrystals()/2);
      if(disk_0) {
        clusterHistD0_->Fill(pos.x(), pos.y(), step_i->energyDepBirks());
      } else {
        clusterHistD1_->Fill(pos.x(), pos.y(), step_i->energyDepBirks());
      }
    }

    if(diagLevel_ > 2
       || (clusterHistD0_->GetMaximum() > minClusterEnergy_ && diagLevel_ > 1)
       || (clusterHistD1_->GetMaximum() > minClusterEnergy_ && diagLevel_ > 1)) {
      std::cout << "Time cluster from " << steps[first_index]->time() << " to " << steps[last_index]->time()
                << " with " << last_index - first_index + 1 << " steps"
                << " Disk 0 max = " << clusterHistD0_->GetMaximum()
                << " Disk 1 max = " << clusterHistD1_->GetMaximum()
                << std::endl;
    }

    // Check for deposits in the histograms above the threshold
    if(checkClusterHistogram(clusterHistD0_)) return true;
    if(checkClusterHistogram(clusterHistD1_)) return true;
    return false;
  }

  //--------------------------------------------------------------------------------
  bool CaloDtsClusterFilter::findCluster(const std::vector<const CaloShowerStep*>& steps) {
    // Loop over steps, looking for time clusters
    const double mbtime = GlobalConstantsHandle<PhysicsParams>()->getNominalDRPeriod();
    size_t first_index(0); // index start for the time cluster
    double edep_cluster(0.); // energy deposit in the cluster
    const size_t nsteps = steps.size();
    for(size_t index = 0; index < nsteps; ++index) {
      const auto* step_i = steps[index];
      edep_cluster += step_i->energyDepBirks();
      if(diagLevel_ > 3) {
        std::cout << "Step " << index
                  << " raw time = " << step_i->time()
                  << " time = " << std::fmod(step_i->time(), mbtime)
                  << " edep = " << step_i->energyDepBirks()
                  << " pdg = " << step_i->simParticle()->pdgId()
                  << std::endl;
      }

      // Check if the first hit is within the cluster, if not move the start to the next step
      const auto* step_first = steps[first_index];
      while(std::fmod(step_i->time(), mbtime) - std::fmod(step_first->time(), mbtime) > timeWindow_) {
        if(first_index == index) break;
        ++first_index;
        edep_cluster -= step_first->energyDepBirks();
        step_first = steps[first_index];
      }

      // Check if the next step is within the cluster, if so continue adding to the cluster
      if(index < nsteps - 1 && (std::fmod(steps[index+1]->time(), mbtime) - std::fmod(step_first->time(), mbtime)) < timeWindow_) {
        continue;
      }

      if(diagLevel_ > 2 || (edep_cluster > minClusterEnergy_ && diagLevel_ > 1)) {
        std::cout << "Step " << index
                  << " time = " << std::fmod(step_i->time(), mbtime)
                  << " edep = " << step_i->energyDepBirks()
                  << " first index = " << first_index
                  << " first time = " << std::fmod(steps[first_index]->time(), mbtime)
                  << " cluster edep = " << edep_cluster
                  << std::endl;
      }

      // Check if the cluster passes the minimum energy cut
      if(edep_cluster < minClusterEnergy_) continue;

      // Next check if the cluster can be clustered in space
      bool in_space_cluster = findSpacialCluster(steps, first_index, index);
      if(in_space_cluster) return true; // found a cluster that passes all cuts, accept
    }
    return false; // no cluster found that passes all cuts, reject
  }

  //--------------------------------------------------------------------------------
  bool CaloDtsClusterFilter::filter(art::Event& event) {
    ++nEvents_;
    if(nullFilter_) {
      ++nPassed_;
      return true;
    }

    // Get the input data
    std::vector<const CaloShowerStep*> steps;
    steps.reserve(1000);
    fillStepsList(steps, event);
    if(steps.empty()) {
      if(diagLevel_ > 1) std::cout << "[CaloDtsClusterFilter::"
        << moduleDescription().moduleLabel()
        << "::" << __func__ << "] " << event.id()
        << " No CaloShowerSteps found\n";
      return false; // no steps, reject
    }

    // Sort steps by time
    sortStepsByTime(steps);

    // Search for clusters within the steps list
    bool result = findCluster(steps);

    // Return the result and update the counter
    if(result) ++nPassed_;
    if(diagLevel_ > 0 && result) {
      std::cout << "[CaloDtsClusterFilter::"
                << moduleDescription().moduleLabel()
                << "::" << __func__ << "] "
                << event.id()
                << " Found cluster passing cuts\n";
    }
    return result;
  }

  //--------------------------------------------------------------------------------
  bool CaloDtsClusterFilter::beginRun(art::Run& run) {
    // Get the calorimeter geometry
    GeomHandle<Calorimeter> ch;
    calorimeter_ = ch.get();
    return true;
  }

  //--------------------------------------------------------------------------------
  void CaloDtsClusterFilter::endJob() {
    mf::LogInfo("Summary")
      << "[CaloDtsClusterFilter::"
      << moduleDescription().moduleLabel()
      << "::" << __func__ << "]"
      << " Passed "
      << nPassed_ << " / " << nEvents_ << " events\n";
  }

}

DEFINE_ART_MODULE(mu2e::CaloDtsClusterFilter)
