///////////////////////////////////////////////////////////////////////////////
// Calorimeter-driven straight line time peak finding
// Michael MacKenzie (2026)
// Based on CalHelixFinder, P. Murat and G. Pezzullo
///////////////////////////////////////////////////////////////////////////////

// framework
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"

// Offline
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/StoppingTargetGeom/inc/StoppingTarget.hh"

#include "Offline/DataProducts/inc/Helicity.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitIndex.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/RecoDataProducts/inc/TrkFitDirection.hh"

//CLHEP
#include "CLHEP/Units/PhysicalConstants.h"

// C++
#include <iostream>
#include <vector>
#include <set>

using CLHEP::HepVector;
using CLHEP::Hep3Vector;

namespace mu2e {

  class CalLineTimePeakFinder : public art::EDProducer {
  protected:

    //-----------------------------------------------------------------------------
    // Main module configuration parameters
    //-----------------------------------------------------------------------------
    struct Config
    {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<int>             diag_level             {Name("DiagLevel"),                  Comment("Diagnostic output level"),0 };
      fhicl::Atom<art::InputTag>   hit_coll_tag           {Name("ComboHitCollectionLabel"),    Comment(" Combo Hit Collection Label") };
      fhicl::Atom<std::string>     calo_cluster_coll_tag  {Name("CaloClusterCollectionLabel"), Comment("CaloCluster Collection Label") };
      fhicl::Atom<int>             min_tc_hits            {Name("MinTimeClusterHits"),         Comment("Min NHits in TimeCluster") };
      fhicl::Atom<float>           min_calo_cluster_energy{Name("MinCaloClusterEnergy"),       Comment("Min Calo Cluster Energy") };
      fhicl::Atom<float>           hit_time_sigma_thresh  {Name("HitTimeSigmaThresh"),         Comment("Time consistency threshold for hits to be added to the cluster (in sigma)") };
      fhicl::Atom<float>           hit_xy_sigma_thresh    {Name("HitXYSigmaThresh"),           Comment("Spatial consistency threshold for hits to be added to the cluster (in sigma)")};
      fhicl::Atom<float>           stopping_target_radius {Name("StoppingTargetRadius"),       Comment("Radius of the stopping target in cone-making (in mm)")};
      fhicl::Atom<float>           particle_beta          {Name("ParticleBeta"),               Comment("Particle beta for distance to time (0-1)"), 1.};
      fhicl::Atom<std::string>     fit_direction          {Name("FitDirection"),               Comment("Fit Direction in Search (\"downstream\" or \"upstream\")") };
    };

    //-----------------------------------------------------------------------------
    // Inputs
    //-----------------------------------------------------------------------------
    art::InputTag                         hit_tag_;                 // input hit collection label
    art::InputTag                         calo_cluster_tag_;        // input calo cluster collection label
    int                                   min_tc_hits_;             // N(hits) in the time cluster
    float                                 min_calo_cluster_energy_; // min energy of the associated calo cluster
    float                                 hit_time_sigma_thresh_;   // time consistency threshold for hits to be added to the cluster
    float                                 hit_xy_sigma_thresh_;     // spatial consistency threshold for hits to be added to the cluster
    float                                 stopping_target_radius_;  // radius of the stopping target in cone-making (in mm)
    float                                 particle_beta_;           // particle beta (v/c)
    TrkFitDirection                       fit_dir_;                 // fit direction in search
    int                                   diag_level_;              // diagnostic output

    //-----------------------------------------------------------------------------
    // Data
    //-----------------------------------------------------------------------------
    const ComboHitCollection*             combo_hit_col_;           // input combo hit collection
    const CaloClusterCollection*          calo_cluster_col_;        // input calo cluster collection
    art::Handle<ComboHitCollection>       combo_hit_col_handle_;    // handle for input combo hit collection
    art::Handle<CaloClusterCollection>    calo_cluster_col_handle_; // handle for input calo cluster collection

    const Tracker*                        tracker_     ;            // tracker geometry
    const Calorimeter*                    calorimeter_ ;            // calorimeter geometry
    const StoppingTarget*                 stopping_target_;         // stopping target geometry
    CLHEP::Hep3Vector                     target_pos_;              // position of the target center for seeding
    double                                calo_d0_offset_;          // z offset of the calorimeter disk 0 from the tracker system
    double                                calo_d1_offset_;          // z offset of the calorimeter disk 1 from the tracker system
  public:
    explicit CalLineTimePeakFinder(const art::EDProducer::Table<Config>& config);
    virtual ~CalLineTimePeakFinder();

    virtual void beginJob();
    virtual void beginRun(art::Run&   run   );
    virtual void produce (art::Event& event );
    virtual void endJob();
    //-----------------------------------------------------------------------------
    // helper functions
    //-----------------------------------------------------------------------------
    bool findData             (const art::Event& e);
    bool isGoodHit            (const ComboHit& hit);
    bool isGoodTimeCluster    (const TimeCluster& tc);
    void findTimePeakInCluster(TimeCluster& tc, const CaloCluster& cl);
    void finalizeTimeCluster  (TimeCluster& tc, const size_t cl_index);
  };

  //-----------------------------------------------------------------------------
  // module constructor, parameter defaults are defined in CalPatRec/fcl/prolog.fcl
  //-----------------------------------------------------------------------------
  CalLineTimePeakFinder::CalLineTimePeakFinder(const art::EDProducer::Table<Config>& config) :
    art::EDProducer{config}
    , hit_tag_(config().hit_coll_tag())
    , calo_cluster_tag_(config().calo_cluster_coll_tag())
    , min_tc_hits_(config().min_tc_hits())
    , min_calo_cluster_energy_(config().min_calo_cluster_energy())
    , hit_time_sigma_thresh_(config().hit_time_sigma_thresh())
    , hit_xy_sigma_thresh_(config().hit_xy_sigma_thresh())
    , stopping_target_radius_(config().stopping_target_radius())
    , particle_beta_(config().particle_beta())
    , fit_dir_(config().fit_direction())
    , diag_level_(config().diag_level())
   {
    // declare the data products
    consumes<ComboHitCollection>(hit_tag_);
    consumes<CaloClusterCollection>(calo_cluster_tag_);
    produces<TimeClusterCollection>();
  }

  //-----------------------------------------------------------------------------
  // destructor
  //-----------------------------------------------------------------------------
  CalLineTimePeakFinder::~CalLineTimePeakFinder() {
  }

  //-----------------------------------------------------------------------------
  void CalLineTimePeakFinder::beginJob() {
  }

  //-----------------------------------------------------------------------------
  void CalLineTimePeakFinder::beginRun(art::Run& ) {
    mu2e::GeomHandle<mu2e::Tracker> th;
    tracker_ = th.get();

    mu2e::GeomHandle<mu2e::Calorimeter> ch;
    calorimeter_ = ch.get();

    mu2e::GeomHandle<mu2e::StoppingTarget> sth;
    stopping_target_ = sth.get();

    // get the offset between the calo disks and the tracker system
    calo_d0_offset_ = calorimeter_->geomUtil().mu2eToTracker(calorimeter_->geomUtil().diskToMu2e(0, CLHEP::Hep3Vector(0., 0., 0.))).z();
    calo_d1_offset_ = calorimeter_->geomUtil().mu2eToTracker(calorimeter_->geomUtil().diskToMu2e(1, CLHEP::Hep3Vector(0., 0., 0.))).z();

    // get the target position in the tracker system
    target_pos_ = calorimeter_->geomUtil().mu2eToTracker(stopping_target_->centerInMu2e());

  }

  //-----------------------------------------------------------------------------
  // find the input data objects
  //-----------------------------------------------------------------------------
  bool CalLineTimePeakFinder::findData(const art::Event& evt) {
    combo_hit_col_ = nullptr;
    calo_cluster_col_ = nullptr;

    evt.getByLabel(hit_tag_, combo_hit_col_handle_);
    if(!combo_hit_col_handle_.isValid()) {
      printf("[CalLineTimePeakFinder::%s] ERROR: ComboHit collection with label \"%s\" not found! RETURN\n", __func__, hit_tag_.encode().c_str());
      combo_hit_col_ = nullptr;
      return false;
    } else {
      combo_hit_col_ = combo_hit_col_handle_.product();
    }

    evt.getByLabel(calo_cluster_tag_, calo_cluster_col_handle_);
    if(!calo_cluster_col_handle_.isValid()) {
      printf("[CalLineTimePeakFinder::%s] ERROR: CaloCluster collection with label \"%s\" not found! RETURN\n", __func__, calo_cluster_tag_.encode().c_str());
      calo_cluster_col_ = nullptr;
      return false;
    } else {
      calo_cluster_col_ = calo_cluster_col_handle_.product();
    }

    return true;
  }

  //-----------------------------------------------------------------------------
  // Find time cluster around a calorimeter cluster
  //-----------------------------------------------------------------------------
  void CalLineTimePeakFinder::findTimePeakInCluster(TimeCluster& tc, const CaloCluster& cl) {

    const float cl_time = cl.time();
    const bool downstream = fit_dir_ == TrkFitDirection::FitDirection::downstream;

    CLHEP::Hep3Vector cl_pos = cl.cog3Vector();
    cl_pos.setZ(cl_pos.z() + (cl.diskID() == 0 ? calo_d0_offset_ : calo_d1_offset_)); // shift the cluster position to the tracker system
    const double dz_target = cl_pos.z() - target_pos_.z(); // z distance from the target to the cluster

    CLHEP::Hep3Vector seed_dir = (downstream) ? cl_pos - target_pos_ : target_pos_ - cl_pos;
    const double seed_dir_mag = seed_dir.mag();
    if(seed_dir_mag <= 0.) return; // can't define a seed direction, so return false{
    seed_dir *= 1./seed_dir_mag; // normalize the seed direction
    if(diag_level_ > 1) {
      printf("[CalLineTimePeakFinder::%s] CaloCluster time = %.2f, position = (%.1f, %.1f, %.1f), seed direction = (%.2f, %.2f, %.2f)\n",
             __func__, cl_time, cl_pos.x(), cl_pos.y(), cl_pos.z(), seed_dir.x(), seed_dir.y(), seed_dir.z());
    }

    // Look for hits consistent this seed position and direction, and add them to the seed
    const size_t n_hits = combo_hit_col_->size();
    for(size_t i_hit = 0; i_hit < n_hits; ++i_hit) {
      const auto& hit = combo_hit_col_->at(i_hit);
      if(!isGoodHit(hit)) continue; // skip hits that don't pass selection
      const CLHEP::Hep3Vector hit_pos(hit.pos().x(), hit.pos().y(), hit.pos().z());
      const double hit_time = hit.correctedTime();
      const double dz = hit_pos.z() - cl_pos.z(); // distance along z to the hit
      const double dr = dz/seed_dir.z(); // distance along the trajectory to the hit location in z
      const double dt = dr / (particle_beta_*CLHEP::c_light); // time of flight to the hit z location
      const double time_at_hit = cl_time + dt; // expected time at the hit z position
      const double time_unc = std::sqrt(hit.timeRes()*hit.timeRes() + cl.timeErr()*cl.timeErr()); // uncertainty on the time difference between the hit and the expected time based on the seed
      const double time_sigma = std::abs((hit_time - time_at_hit) / time_unc); // number of sigma the hit time is from the expected time based on the seed
      if(diag_level_ > 2) {
        printf(" Hit %zu: time = %.2f, expected time = %.2f, time sigma = %.2f\n",
               i_hit, hit_time, time_at_hit, time_sigma);
      }
      if(time_sigma > hit_time_sigma_thresh_) continue; // hit is not consistent with the seed, so skip it

      // next check if the hit is consistent in space
      const CLHEP::Hep3Vector pos_at_hit = cl_pos + dr * seed_dir; // position along the seed direction at the same z as the hit
      const double x_y_dist = (hit_pos - pos_at_hit).perp(); // distance in the x-y plane between the hit and the expected position based on the seed
      const double cone = std::abs(dz/dz_target) * stopping_target_radius_; // radius of the cone in the x-y plane at the hit z based on the seed direction and cluster position
      const double xy_unc = std::sqrt(hit.transVar() + hit.wireVar()); // uncertainty in the hit position in the x-y plane
      const double xy_sigma = (x_y_dist - cone) / xy_unc; // number of sigma the hit is from the cone from the cluster to the target, negative means contained within the cone
      if(diag_level_ > 2) {
        printf("  Expected position at hit z: (%.1f, %.1f, %.1f), hit position: (%.1f, %.1f, %.1f), x-y distance = %.2f, target cone = %.2f, sigma = %.2f\n",
               pos_at_hit.x(), pos_at_hit.y(), pos_at_hit.z(), hit_pos.x(), hit_pos.y(), hit_pos.z(), x_y_dist, cone, xy_sigma);
      }
      if(xy_sigma > hit_xy_sigma_thresh_) continue; // hit is not consistent with the seed
      tc._strawHitIdxs.push_back(i_hit);
    }
  }

  //-----------------------------------------------------------------------------
  // event entry point
  //-----------------------------------------------------------------------------
  void CalLineTimePeakFinder::produce(art::Event& event ) {

    // output collection
    std::unique_ptr<TimeClusterCollection> out_tcs(new TimeClusterCollection);

    //-----------------------------------------------------------------------------
    // find the data
    //-----------------------------------------------------------------------------

    const bool valid_data = findData(event);

    //-----------------------------------------------------------------------------
    // Search around each calo cluster for time peaks
    //-----------------------------------------------------------------------------

    if(valid_data) {
      const size_t n_clusters = calo_cluster_col_->size();
      for(size_t i_cl = 0; i_cl < n_clusters; ++i_cl) {
        const auto& cl = calo_cluster_col_->at(i_cl);
        if(diag_level_ > 1) {
          printf("[CalLineTimePeakFinder::%s] Seeding from CaloCluster %zu with energy = %.2f, time = %.2f, N(combo hits) = %zu\n",
                 __func__, i_cl, cl.energyDep(), cl.time(), combo_hit_col_->size());
        }
        if(cl.energyDep() < min_calo_cluster_energy_) continue; // skip clusters that don't pass the energy cut
        TimeCluster tc;
        findTimePeakInCluster(tc, cl);
        if(!isGoodTimeCluster(tc)) continue; // skip time clusters that don't pass selection criteria
        finalizeTimeCluster(tc, i_cl);
        out_tcs->push_back(std::move(tc));
      }
    }

    //-----------------------------------------------------------------------------
    // put reconstructed time clusters into the event record
    //-----------------------------------------------------------------------------

    if(diag_level_ > 0) {
      printf("[CalLineTimePeakFinder::%s] Found %zu time clusters in total\n", __func__, out_tcs->size());
    }
    event.put(std::move(out_tcs));
  }

  //-----------------------------------------------------------------------------
  bool  CalLineTimePeakFinder::isGoodHit(const ComboHit& hit) {
    const auto&  flag    = hit.flag();
    const int    bkg_hit = flag.hasAnyProperty(StrawHitFlag::bkg);
    if (bkg_hit != 0) return false;
    return true;
  }

  //-----------------------------------------------------------------------------
  bool CalLineTimePeakFinder::isGoodTimeCluster(const TimeCluster& tc) {
    if(int(tc.nhits()) < min_tc_hits_) return false; // not enough hits in the time cluster
    return true;
  }

  //-----------------------------------------------------------------------------
  void CalLineTimePeakFinder::finalizeTimeCluster(TimeCluster& tc, size_t cl_index) {
    tc._caloCluster = art::Ptr<CaloCluster>(calo_cluster_col_handle_, cl_index);
    XYZVectorF tc_pos(0.,0.,0.);
    double avg_t0 = 0.;
    double t0sq = 0.;
    double weight = 0.;
    for(size_t hit_index : tc._strawHitIdxs) {
      const auto& hit = combo_hit_col_->at(hit_index);
      tc._nsh += hit.nStrawHits();
      tc_pos += hit.pos() * hit.nStrawHits();
      avg_t0 += hit.correctedTime() * hit.nStrawHits();
      t0sq += std::pow(hit.correctedTime(), 2) * hit.nStrawHits();
      weight += hit.nStrawHits();
    }
    if(weight > 0.) {
      avg_t0 /= weight;
      t0sq /= weight;
      tc_pos /= weight;
    }
    tc._pos = tc_pos;
    tc._t0._t0 = avg_t0;
    tc._t0._t0err = std::sqrt(std::max(0., t0sq  - avg_t0*avg_t0));
  }

  //-----------------------------------------------------------------------------
  //
  //-----------------------------------------------------------------------------
  void CalLineTimePeakFinder::endJob() {
  }

} // end namespace mu2e

using mu2e::CalLineTimePeakFinder;
DEFINE_ART_MODULE(CalLineTimePeakFinder)
