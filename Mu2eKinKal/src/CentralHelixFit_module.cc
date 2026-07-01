//
// KinKal fit module using an input helix seed.  Base for either central or loop helix fit
//
// Original author D. Brown (LBNL) 11/18/20
//
// framework
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Tuple.h"
#include "fhiclcpp/types/OptionalTable.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
// conditions
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/TrackerConditions/inc/StrawResponse.hh"
#include "Offline/BFieldGeom/inc/BFieldManager.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
// utiliites
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeneralUtilities/inc/Angles.hh"
#include "Offline/TrkReco/inc/TrkUtilities.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/GeneralUtilities/inc/OwningPointerCollection.hh"
#include "Offline/Mu2eUtilities/inc/TwoLinePCA_XYZ.hh"
// data
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/DataProducts/inc/Helicity.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"
#include "Offline/RecoDataProducts/inc/KalSeedAssns.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/CosmicTrackSeed.hh"
#include "Offline/DataProducts/inc/SurfaceId.hh"
#include "Offline/KinKalGeom/inc/KinKalGeom.hh"
// MC (TruthSeedDiag only)
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/MCTrajectory.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
// KinKal
#include "KinKal/Fit/Track.hh"
#include "KinKal/Fit/Config.hh"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "KinKal/Trajectory/PiecewiseClosestApproach.hh"
#include "KinKal/Trajectory/CentralHelix.hh"
#include "KinKal/General/Parameters.hh"
#include "KinKal/General/Vectors.hh"
// Mu2eKinKal
#include "Offline/Mu2eKinKal/inc/KKFit.hh"
#include "Offline/Mu2eKinKal/inc/KKFitSettings.hh"
#include "Offline/Mu2eKinKal/inc/KKTrack.hh"
#include "Offline/Mu2eKinKal/inc/KKStrawHit.hh"
#include "Offline/Mu2eKinKal/inc/KKStrawHitCluster.hh"
#include "Offline/Mu2eKinKal/inc/KKStrawXing.hh"
#include "Offline/Mu2eKinKal/inc/KKCaloHit.hh"
#include "Offline/Mu2eKinKal/inc/KKBField.hh"
#include "Offline/Mu2eKinKal/inc/KKConstantBField.hh"
#include "Offline/Mu2eKinKal/inc/KKFitUtilities.hh"
#include "Offline/Mu2eKinKal/inc/KKExtrap.hh"
// root
#include "TH1F.h"
#include "TTree.h"
#include "Math/AxisAngle.h"
// C++
#include <iostream>
#include <fstream>
#include <string>
#include <functional>
#include <vector>
#include <memory>
#include <array>
#include <cmath>
#include <limits>

using KTRAJ= KinKal::CentralHelix; // this must come before HelixFit
using namespace ROOT::Math;
#include "Offline/TrkReco/inc/TrkUtilities.hh"
#include "Offline/GeneralUtilities/inc/Angles.hh"

namespace mu2e {
  using PTRAJ = KinKal::ParticleTrajectory<KTRAJ>;
  using KKTRK = KKTrack<KTRAJ>;
  using KKTRKCOL = OwningPointerCollection<KKTRK>;
  using KKSTRAWHIT = KKStrawHit<KTRAJ>;
  using KKSTRAWHITPTR = std::shared_ptr<KKSTRAWHIT>;
  using KKSTRAWHITCOL = std::vector<KKSTRAWHITPTR>;
  using KKSTRAWXING = KKStrawXing<KTRAJ>;
  using KKSTRAWXINGPTR = std::shared_ptr<KKSTRAWXING>;
  using KKSTRAWXINGCOL = std::vector<KKSTRAWXINGPTR>;
  using KKCALOHIT = KKCaloHit<KTRAJ>;
  using KKCALOHITPTR = std::shared_ptr<KKCALOHIT>;
  using KKCALOHITCOL = std::vector<KKCALOHITPTR>;
  using PARAMHIT = KinKal::ParameterHit<KTRAJ>;
  using PARAMHITPTR = std::shared_ptr<PARAMHIT>;
  using PARAMHITCOL = std::vector<PARAMHITPTR>;
  using KKFIT = KKFit<KTRAJ>;
  using KinKal::VEC3;
  using KinKal::DMAT;
  using HPtr = art::Ptr<HelixSeed>;
  using CCPtr = art::Ptr<CaloCluster>;
  using CCHandle = art::Handle<CaloClusterCollection>;
  using StrawHitIndexCollection = std::vector<StrawHitIndex>;

  using KKConfig = Mu2eKinKal::KinKalConfig;
  using KKFitConfig = Mu2eKinKal::KKFitConfig;
  using KKModuleConfig = Mu2eKinKal::KKModuleConfig;

  using Name    = fhicl::Name;
  using Comment = fhicl::Comment;
  using KinKal::DVEC;

  // extend the generic module configuration as needed
  struct KKCHModuleConfig : KKModuleConfig {
    fhicl::Sequence<art::InputTag> seedCollections         {Name("CosmicTrackSeedCollections"),     Comment("Seed fit collections to be processed ") };
    fhicl::OptionalAtom<double> fixedBField { Name("ConstantBField"), Comment("Constant BField value") };
    fhicl::Atom<double> seedMom { Name("SeedMomentum"), Comment("Seed momentum") };
    fhicl::Atom<int> seedCharge { Name("SeedCharge"), Comment("Seed charge MAGNITUDE, in electron charge units") };
    fhicl::Sequence<double> parconst { Name("ParameterConstraints"), Comment("External constraint on parameters to seed values (rms, various units)") };
    fhicl::Sequence<std::string> sampleSurfaces { Name("SampleSurfaces"), Comment("When creating the KalSeed, sample the fit at these surfaces") };
    fhicl::Atom<bool> sampleInRange { Name("SampleInRange"), Comment("Require sample times to be inside the fit trajectory time range") };
    fhicl::Atom<bool> sampleInBounds { Name("SampleInBounds"), Comment("Require sample intersection point be inside surface bounds (within tolerance)") };
    fhicl::Atom<float> interTol { Name("IntersectionTolerance"), Comment("Tolerance for surface intersections (mm)") };
    fhicl::Atom<float> sampleTBuff { Name("SampleTimeBuffer"), Comment("Time buffer for sample intersections (nsec)") };
    fhicl::Atom<bool> useFitCharge { Name("UseFitCharge"), Comment("Set the PDG particle according to the fit charge; otherwise reject fits that don't agree with the PDG particle charge") };
    fhicl::Atom<float> minCenterRho { Name("MinCenterRho"), Comment("Minimum transverse distance from the helix axis to the Z axis to consider the fit non-degenerate (mm)") };
    fhicl::Atom<bool> truthSeedDiag { Name("TruthSeedDiag"), Comment("Diagnostic: also extrapolate a copy of each fit re-seeded with the true muon state, isolating extrapolation error from fit error. Writes KalSeeds under instance 'TruthSeed'"), false };
    fhicl::OptionalAtom<art::InputTag> mcTrajectoryCollection { Name("MCTrajectories"), Comment("MCTrajectory map providing the true muon state at the fit reference point (required iff TruthSeedDiag)") };
    fhicl::Sequence<double> truthSeedParamConstraints {
      Name("TruthSeedParameterConstraints"),
      Comment("TruthSeedDiag ParameterHit RMS constraints for CentralHelix parameters d0, phi0, omega, z0, tanDip, t0"),
      std::vector<double>{1.0e-3, 1.0e-6, 1.0e-9, 1.0e-3, 1.0e-6, 1.0e-3}
    };
  };

  struct GlobalConfig {
    fhicl::Table<KKCHModuleConfig> modSettings { Name("ModuleSettings") };
    fhicl::Table<KKFitConfig> kkfitSettings { Name("KKFitSettings") };
    fhicl::Table<KKConfig> fitSettings { Name("FitSettings") };
    fhicl::Table<KKConfig> extSettings { Name("ExtensionSettings") };
    fhicl::OptionalTable<KKExtrapConfig> extrapSettings { Name("ExtrapolationSettings") };
    // helix module specific config
  };

  class CentralHelixFit : public art::EDProducer {
    public:
      using Parameters = art::EDProducer::Table<GlobalConfig>;
      explicit CentralHelixFit(const Parameters& settings);
      virtual ~CentralHelixFit() {}
      void beginRun(art::Run& run) override;
      void produce(art::Event& event) override;
    protected:
      bool goodFit(KKTRK const& ktrk) const;
      void sampleFit(KKTRK& ktrk) const;
      TrkFitFlag fitflag_;
      // parameter-specific functions that need to be overridden in subclasses
      // data payload
      art::ProductToken<ComboHitCollection> chcol_T_;
      art::ProductToken<CaloClusterCollection> cccol_T_;
      std::vector<art::ProductToken<CosmicTrackSeedCollection>> cseedCols_;
      TrkFitFlag goodseed_;
      bool saveall_;
      ProditionsHandle<StrawResponse> strawResponse_h_;
      ProditionsHandle<Tracker> alignedTracker_h_;
      int print_;
      PDGCode::type fpart_;
      KKFIT kkfit_; // fit helper
      DMAT seedcov_; // seed covariance matrix
      double mass_; // particle mass
      int PDGcharge_; // PDG particle charge
      std::unique_ptr<KinKal::BFieldMap> kkbf_;
      Config config_; // initial fit configuration object
      Config exconfig_; // extension configuration object
      std::unique_ptr<KKExtrap> extrap_; // extrapolation helper
      bool fixedfield_; //
      double seedMom_;
      int seedCharge_;
      std::vector<double> paramconstraints_;
      double intertol_; // surface intersection tolerance (mm)
      double sampletbuff_; // simple time buffer; replace this with extrapolation TODO
      bool useFitCharge_; // Set the PDG particle to agree with the fit charge
      double minCenterRho_; // min center distance to z axis
      bool sampleinrange_, sampleinbounds_; // require samples to be in range or on surface
      SurfaceIdCollection ssids_;
      KinKalGeom::SurfacePairCollection surfacess_to_sample_; // surfaces to sample the fit
      bool constraining_;
      bool truthSeedDiag_; // diagnostic: extrapolate a truth-seeded copy of each fit
      art::ProductToken<MCTrajectoryCollection> mctrajcol_T_; // MC trajectories (truth seed)
      std::vector<double> truthSeedParamConstraints_; // ParameterHit RMS constraints for truth seeding
  };

  CentralHelixFit::CentralHelixFit(const Parameters& settings) : art::EDProducer{settings},
    fitflag_(TrkFitFlag::KKCentralHelix),
    chcol_T_(consumes<ComboHitCollection>(settings().modSettings().comboHitCollection())),
    cccol_T_(mayConsume<CaloClusterCollection>(settings().modSettings().caloClusterCollection())),
    goodseed_(settings().modSettings().seedFlags()),
    saveall_(settings().modSettings().saveAll()),
    print_(settings().modSettings().printLevel()),
    fpart_(static_cast<PDGCode::type>(settings().modSettings().fitParticle())),
    kkfit_(settings().kkfitSettings()),
    config_(Mu2eKinKal::makeConfig(settings().fitSettings())),
    exconfig_(Mu2eKinKal::makeConfig(settings().extSettings())),
    fixedfield_(false),
    seedMom_(settings().modSettings().seedMom()),
    seedCharge_(settings().modSettings().seedCharge()),
    paramconstraints_(settings().modSettings().parconst()),
    intertol_(settings().modSettings().interTol()),
    sampletbuff_(settings().modSettings().sampleTBuff()),
    useFitCharge_(settings().modSettings().useFitCharge()),
    minCenterRho_(settings().modSettings().minCenterRho()),
    sampleinrange_(settings().modSettings().sampleInRange()),
    sampleinbounds_(settings().modSettings().sampleInBounds()),
    truthSeedDiag_(settings().modSettings().truthSeedDiag()),
    truthSeedParamConstraints_(settings().modSettings().truthSeedParamConstraints())
    {
      // collection handling
      for(const auto& cseedtag : settings().modSettings().seedCollections()) { cseedCols_.emplace_back(consumes<CosmicTrackSeedCollection>(cseedtag)); }
      produces<KKTRKCOL>();
      produces<KalSeedCollection>();
      //      produces<KalHelixAssns>();
      // TruthSeedDiag: a parallel truth-seeded extrapolation, written under instance "TruthSeed"
      if(truthSeedDiag_){
        art::InputTag mctrajtag;
        if(!settings().modSettings().mcTrajectoryCollection(mctrajtag))
          throw cet::exception("RECO") << "mu2e::CentralHelixFit: TruthSeedDiag requires MCTrajectories" << endl;
        if(truthSeedParamConstraints_.size() != KinKal::NParams())
          throw cet::exception("RECO") << "mu2e::CentralHelixFit: TruthSeedParameterConstraints must have "
            << KinKal::NParams() << " entries" << endl;
        for(auto const& sigma : truthSeedParamConstraints_){
          if(sigma <= 0.0)
            throw cet::exception("RECO") << "mu2e::CentralHelixFit: TruthSeedParameterConstraints entries must be positive" << endl;
        }
        mctrajcol_T_ = consumes<MCTrajectoryCollection>(mctrajtag);
        produces<KKTRKCOL>("TruthSeed");
        produces<KalSeedCollection>("TruthSeed");
      }
      // build the initial seed covariance
      auto const& seederrors = settings().modSettings().seederrors();
      if(seederrors.size() != KinKal::NParams())
        throw cet::exception("RECO")<<"mu2e::CentralHelixFit:Seed error configuration error"<< endl;
      for(size_t ipar=0;ipar < seederrors.size(); ++ipar){
        seedcov_[ipar][ipar] = seederrors[ipar]*seederrors[ipar];
      }
      constraining_ = false;
      if (paramconstraints_.size() == KinKal::NParams()){
        for (size_t ipar=0;ipar<KinKal::NParams();ipar++){
          if (paramconstraints_[ipar] > 0)
            constraining_ = true;
        }
      }else if (paramconstraints_.size() > 0){
        throw cet::exception("RECO")<<"mu2e::CentralHelixFit: Parameter constraint configuration error"<< endl;
      }

      if(print_ > 0) std::cout << "Fit " << config_ << "Extension " << exconfig_;
      double bz(0.0);
      if(settings().modSettings().fixedBField(bz)){
        fixedfield_ = true;
        kkbf_ = std::move(std::make_unique<KKConstantBField>(VEC3(0.0,0.0,bz)));
      }
      // setup extrapolation
      if(settings().extrapSettings())extrap_ = make_unique<KKExtrap>(*settings().extrapSettings());

      // surfaces to sample; this interface is deprecatecd and should be replaced with extrapolation TODO
      for(auto const& sidname : settings().modSettings().sampleSurfaces()) {
        ssids_.push_back(SurfaceId(sidname,-1)); // match all elements
      }
    }

  void CentralHelixFit::beginRun(art::Run& run) {
    // setup things that rely on data related to beginRun
    auto const& ptable = GlobalConstantsHandle<ParticleDataList>();
    mass_ = ptable->particle(fpart_).mass();
    PDGcharge_ = static_cast<int>(ptable->particle(fpart_).charge());
    // create KKBField
    if(!fixedfield_){
      GeomHandle<BFieldManager> bfmgr;
      GeomHandle<DetectorSystem> det;
      kkbf_ = std::move(std::make_unique<KKBField>(*bfmgr,*det));
    }
    if(print_ > 0) kkbf_->print(std::cout);
    // translate the sample surface names to actual surfaces using the KinKalGeom. This must be done after construction as the KKGeom object now comes from GeometryService
    GeomHandle<mu2e::KinKalGeom> kkg_h;
    auto const& kkg = *kkg_h;
    kkg.surfaces(ssids_,surfacess_to_sample_);
  }

  void CentralHelixFit::produce(art::Event& event ) {
    GeomHandle<Calorimeter> calo_h;
    GeomHandle<mu2e::Tracker> nominalTracker_h;
    // find current proditions
    auto const& strawresponse = strawResponse_h_.getPtr(event.id());
    auto const& tracker = alignedTracker_h_.getPtr(event.id()).get();
    // find input hits
    auto ch_H = event.getValidHandle<ComboHitCollection>(chcol_T_);
    auto cc_H = event.getHandle<CaloClusterCollection>(cccol_T_);
    auto const& chcol = *ch_H;
    // create output
    unique_ptr<KKTRKCOL> ktrkcol(new KKTRKCOL );
    unique_ptr<KalSeedCollection> kkseedcol(new KalSeedCollection );
    // TruthSeedDiag parallel outputs + truth inputs
    unique_ptr<KKTRKCOL> ktrkcol_truth(new KKTRKCOL );
    unique_ptr<KalSeedCollection> kkseedcol_truth(new KalSeedCollection );
    MCTrajectoryCollection const* mctrajs = nullptr;
    GeomHandle<DetectorSystem> det;
    if(truthSeedDiag_) mctrajs = event.getValidHandle<MCTrajectoryCollection>(mctrajcol_T_).product();

    unsigned nseed(0);
    for (auto const& cseedtag : cseedCols_) {
      auto const& cseedcol_h = event.getValidHandle<CosmicTrackSeedCollection>(cseedtag);
      auto const& cseedcol = *cseedcol_h;
      nseed += cseedcol.size();

      for (size_t index=0;index<cseedcol.size();++index){
        const auto& cseed = cseedcol[index];

        auto trange = Mu2eKinKal::timeBounds(cseed.hits());
        auto fitpart = fpart_;


        XYZVectorF trackmid(0,0,0);
        for (size_t i=0;i<cseed.hits().size();i++){
          auto hiti = cseed.hits()[i];
          trackmid += hiti.pos();
        }
        trackmid /= cseed.hits().size();

        XYZVectorF trackpos = cseed.track().FitEquation.Pos;
        XYZVectorF trackdir = cseed.track().FitEquation.Dir.Unit();
        if (trackdir.y() > 0)
          trackdir *= -1;

        // put direction as tangent to circle at mid y position
        double dist = 0;
        if (fabs(trackdir.y()) > 0.01)
          dist = (trackmid.y() - trackpos.y())/trackdir.y();
        else if (fabs(trackdir.x()) > 0.01)
          dist = (trackmid.x() - trackpos.x())/trackdir.x();
        trackpos += dist * trackdir;
        double t0 = cseed.t0().t0() + dist/299.9;

        // take the magnetic field at track center as nominal
        auto bnom = kkbf_->fieldVect(VEC3(trackpos.x(),trackpos.y(),trackpos.z()));

        XYZVectorF trackmom = trackdir * seedMom_;

        KinKal::Parameters kkpars;
        try {
          auto temptraj = KTRAJ(KinKal::VEC4(trackpos.x(),trackpos.y(),trackpos.z(),t0),KinKal::MOM4(trackmom.x(),trackmom.y(),trackmom.z(),mass_), seedCharge_, bnom.Z());
          kkpars = KinKal::Parameters(temptraj.params().parameters(),seedcov_);
        } catch (std::invalid_argument const& error) {
          if(print_ > 0) std::cout << "CentralHelixFit Seed Error " << error.what() << std::endl;
          continue;
        }
        auto seedtraj = KTRAJ(kkpars,mass_,seedCharge_,bnom.Z(),trange);

        // wrap the seed traj in a Piecewise traj: needed to satisfy PTOCA interface
        PTRAJ pseedtraj(seedtraj);

        // first, we need to unwind the combohits.  We use this also to find the time range
        StrawHitIndexCollection strawHitIdxs;
        auto chcolptr = cseed.hits().fillStrawHitIndices(strawHitIdxs, StrawIdMask::uniquestraw);
        if(chcolptr != &chcol)
          throw cet::exception("RECO")<<"mu2e::KKCentralHelixFit: inconsistent ComboHitCollection" << std::endl;
        // next, build straw hits and materials from these
        KKSTRAWHITCOL strawhits;
        strawhits.reserve(strawHitIdxs.size());
        KKSTRAWXINGCOL strawxings;
        strawxings.reserve(strawHitIdxs.size());
        kkfit_.makeStrawHits(*tracker, *strawresponse, *kkbf_, pseedtraj, chcol, strawHitIdxs, strawhits, strawxings);
        // optionally (and if present) add the CaloCluster as a constraint
        // verify the cluster looks physically reasonable before adding it TODO!  Or, let the KKCaloHit updater do it TODO
        KKCALOHITCOL calohits;
        //FIXME    if (kkfit_.useCalo() && hseed.caloCluster().isNonnull())kkfit_.makeCaloHit(hseed.caloCluster(),*calo_h, pseedtraj, calohits);
        // extend the seed range given the hits and xings

        try {
          seedtraj.range() = kkfit_.range(strawhits,calohits,strawxings);
          // create parameter constraint, use seedtrajectory parameters
          PARAMHITCOL paramhits;
          if (constraining_){
            kkfit_.makeSeedParamHit(seedtraj,paramconstraints_,paramhits);
          }

          // create and fit the track
          auto ktrk = make_unique<KKTRK>(config_,*kkbf_,seedtraj,fitpart,kkfit_.strawHitClusterer(),strawhits,strawxings,calohits,paramhits);
          // Check the fit
          auto goodfit = goodFit(*ktrk);
          // if we have an extension schedule, extend.
          if(goodfit && exconfig_.schedule().size() > 0) {
            //  std::cout << "EXTENDING TRACK " << event.id() << " " << index << std::endl;
            kkfit_.extendTrack(exconfig_,*kkbf_, *tracker,*strawresponse, chcol, *calo_h, cc_H, *ktrk );
            goodfit = goodFit(*ktrk);
          }
          // [TruthSeedDiag] Build a parallel track from the piecewise tracker truth trajectory,
          // constrained to it by ParameterHits, and extrapolate that (via the PUBLIC KinKal
          // reconstruction interface -- no fitTrajMutable substitution). The truth-seeded CRV
          // residual isolates the extrapolation (B-field/material) error; differencing it from
          // the fit-seeded residual gives the fit's contribution (the charge over-rotation).
          if(truthSeedDiag_ && goodfit && extrap_ && mctrajs != nullptr){
            // The fit trajectory is PIECEWISE (one CentralHelix per domain) with dE/dx applied
            // between pieces, so its |p| steps down through the tracker. A single truth helix has
            // CONSTANT |p| and so cannot match the truth at every tracker surface for a low-|p|
            // muon that sheds a large momentum fraction crossing the tracker. So we rebuild the
            // truth seed with the SAME domain structure as the fit: one truth helix per fit piece,
            // each anchored at the true muon state (MC trajectory) at that piece's location. This
            // reproduces the truth |p| -- tracker energy loss included -- at every surface, giving a
            // ~zero tracker residual at all momenta. The extrapolation then carries it to the CRV.
            // NB: do NOT dereference the trajectory's SimParticle Ptr -- cosmic resampling/mixing
            // leaves it dangling (ProductNotFound); the trajectory POINTS (pos/time/KE) are valid.
            int tcharge = (ktrk->fitTraj().nearestPiece(ktrk->fitTraj().t0()).charge() > 0.0) ? 1 : -1; // from the fit
            // truth muon state (pos, mom, time) from the MC trajectory point nearest a given location
            auto trueStateAt = [&](VEC3 const& at, VEC3& tpos, VEC3& tmom, double& ttime)->bool{
              double bestd2 = std::numeric_limits<double>::max(); bool found = false;
              for(auto const& tjpair : *mctrajs){
                auto const& pts = tjpair.second.points();
                if(pts.size() < 2) continue;
                for(size_t i = 0; i < pts.size(); ++i){
                  auto pdh = det->toDetector(pts[i].pos());
                  VEC3 pd(pdh.x(),pdh.y(),pdh.z());
                  double d2 = (pd - at).Mag2();
                  if(d2 >= bestd2) continue;
                  size_t j = (i + 1 < pts.size()) ? i + 1 : i - 1; // local tangent for the direction
                  auto qdh = det->toDetector(pts[j].pos());
                  VEC3 qd(qdh.x(),qdh.y(),qdh.z());
                  VEC3 dir = (j > i) ? (qd - pd) : (pd - qd);
                  double dn = sqrt(dir.Mag2());
                  if(dn <= 0.0) continue;
                  dir /= dn;
                  double ke = pts[i].kineticEnergy();
                  double pmag = std::sqrt(ke*ke + 2.0*ke*mass_); // |p| from kinetic energy
                  bestd2 = d2; found = true; tpos = pd; tmom = pmag*dir; ttime = pts[i].t();
                }
              }
              return found;
            };
            // build the piecewise truth trajectory (one truth helix per fit domain) and remember each
            // piece's reference time, to constrain it with a ParameterHit below.
            PTRAJ truthpt;
            std::vector<double> truthHitTimes;
            for(auto const& fpieceptr : ktrk->fitTraj().pieces()){
              auto const& fpiece = *fpieceptr;
              auto const htime = fpiece.range().mid();
              VEC3 tpos, tmom; double ttime = 0.0;
              if(!trueStateAt(fpiece.position3(htime), tpos, tmom, ttime)) continue;
              auto tbnom = kkbf_->fieldVect(tpos);
              KTRAJ tpiece(KinKal::VEC4(tpos.X(),tpos.Y(),tpos.Z(),ttime),
                           KinKal::MOM4(tmom.X(),tmom.Y(),tmom.Z(),mass_),
                           tcharge, tbnom.Z(), fpiece.range());
              tpiece.params() = KinKal::Parameters(tpiece.params().parameters(), seedcov_);
              truthpt.append(tpiece);
              truthHitTimes.emplace_back(htime);
            }
            if(!truthHitTimes.empty()){
              // Constrain each truth piece with a ParameterHit at its reference time (RMS from
              // TruthSeedParameterConstraints). This pins the fit to the truth via the PUBLIC KinKal
              // reconstruction interface -- no fitTrajMutable substitution -- at the cost of being a
              // (tightly) constrained fit rather than an exact trajectory swap.
              PARAMHITCOL truthParamHits;
              PARAMHIT::PMASK truthMask;
              truthMask.fill(true);
              auto addTruthHit = [&](double hitTime, double covarianceScale) {
                auto cparams = truthpt.nearestPiece(hitTime).params();
                for(size_t ipar = 0; ipar < KinKal::NParams(); ++ipar){
                  for(size_t jpar = 0; jpar < KinKal::NParams(); ++jpar){
                    cparams.covariance()[ipar][jpar] = 0.0;
                  }
                  auto const sigma = truthSeedParamConstraints_.at(ipar);
                  cparams.covariance()[ipar][ipar] = covarianceScale*sigma*sigma;
                }
                truthParamHits.push_back(std::make_shared<PARAMHIT>(
                    hitTime, truthpt, cparams, truthMask));
              };
              // a single piece gives a zero-length fit range (lowNDOF); split it into two hits a
              // sliver apart, with looser (x2) constraint, as the line tail does.
              if(truthHitTimes.size() == 1){
                constexpr double truthHitDt = 1.0e-3;
                addTruthHit(truthHitTimes.front() - 0.5*truthHitDt, 2.0);
                addTruthHit(truthHitTimes.front() + 0.5*truthHitDt, 2.0);
              } else {
                for(auto const hitTime : truthHitTimes) addTruthHit(hitTime, 1.0);
              }
              // copy the fit's field domains so the truth fit/extrapolation see the same BField model
              KKTRK::DOMAINCOL truthDomains;
              for(auto const& domain : ktrk->domains()){
                truthDomains.emplace(std::make_shared<KinKal::Domain>(*domain));
              }
              KKTRK::PKTRAJPTR truthTraj = std::make_unique<PTRAJ>(truthpt);
              KKSTRAWHITCOL nostrawhits; KKSTRAWXINGCOL nostrawxings; KKCALOHITCOL nocalohits;
              // one iteration, no convergence churn: the ParameterHits already pin the params
              auto truthConfig = config_;
              truthConfig.maxniter_ = 1;
              truthConfig.divdchisq_ = std::numeric_limits<double>::max();
              truthConfig.pdchisq_ = std::numeric_limits<double>::max();
              truthConfig.divgap_ = std::numeric_limits<double>::max();
              auto ktrk_truth = make_unique<KKTRK>(truthConfig,*kkbf_,fitpart,truthTraj,
                  nostrawhits,nostrawxings,nocalohits,truthParamHits,truthDomains);
              if(ktrk_truth->fitStatus().usable()){
                extrap_->extrapolate(*ktrk_truth);
                TrkFitFlag tflag; tflag.merge(fitflag_); tflag.merge(TrkFitFlag::FitOK);
                auto truthKSeed = kkfit_.createSeed(*ktrk_truth,tflag,*calo_h,*nominalTracker_h);
                // EventNtuple and other consumers use the detector-hit list for bookkeeping even
                // when hit details are disabled; carry the reco fit's hit metadata onto the truth
                // KalSeed without adding the detector measurements to the truth-constrained fit.
                auto recoMetadata = kkfit_.createSeed(*ktrk,tflag,*calo_h,*nominalTracker_h);
                truthKSeed._hits = std::move(recoMetadata._hits);
                truthKSeed._hitcalibs = std::move(recoMetadata._hitcalibs);
                truthKSeed._chit = std::move(recoMetadata._chit);
                kkseedcol_truth->push_back(std::move(truthKSeed));
                ktrkcol_truth->push_back(ktrk_truth.release());
              } else if(print_ > 0) {
                std::cout << "CentralHelixFit TruthSeed ParameterHit fit unusable: "
                  << ktrk_truth->fitStatus() << std::endl;
              }
            }
          }
          // extrapolate as required
          if(goodfit && extrap_)extrap_->extrapolate(*ktrk);
          if(print_>1)ktrk->printFit(std::cout,print_);
          if(goodfit || saveall_){
            TrkFitFlag fitflag;
            fitflag.merge(fitflag_);
            if(goodfit)
              fitflag.merge(TrkFitFlag::FitOK);
            else
              fitflag.clear(TrkFitFlag::FitOK);
            // flip the PDG particle assignment charge if required
            double t0charge = ktrk->fitTraj().nearestPiece(ktrk->fitTraj().t0()).charge();
            if(t0charge*PDGcharge_ < 0)ktrk->reverseCharge();
            // sample as requested: this may be redundant with extrapolation
            sampleFit(*ktrk);
            auto kkseed = kkfit_.createSeed(*ktrk,fitflag,*calo_h,*nominalTracker_h);
            kkseedcol->push_back(kkseed);
            ktrkcol->push_back(ktrk.release());
          }
        } catch (std::invalid_argument const& error) {
          if(print_ > 0) std::cout << "CentralHelixFit Error " << error.what() << std::endl;
        }
      }
    }
    // put the output products into the event
    if(print_ > 0) std::cout << "Fitted " << ktrkcol->size() << " tracks from " << nseed << " Seeds" << std::endl;
    event.put(move(ktrkcol));
    event.put(move(kkseedcol));
    if(truthSeedDiag_){
      event.put(move(ktrkcol_truth),"TruthSeed");
      event.put(move(kkseedcol_truth),"TruthSeed");
    }
  }

  bool CentralHelixFit::goodFit(KKTRK const& ktrk) const {
    // require physical consistency: fit can succeed but the result can have changed charge or helicity
    bool retval = ktrk.fitStatus().usable();
    // compare fit charge with intent
    if(retval){
      double t0charge = ktrk.fitTraj().nearestPiece(ktrk.fitTraj().t0()).charge();
      retval = t0charge*PDGcharge_> 0 || useFitCharge_;
    }
    // also check that the fit is inside the physical detector volume.  Test where the StrawHits are
    if(retval){
      for(auto const& shptr : ktrk.strawHits()) {
        if(shptr->active() && !Mu2eKinKal::inDetector(ktrk.fitTraj().position3(shptr->time()))){
          retval = false;
          break;
        }
      }
    }
    if(retval){
      // check that the segments are non-degenerate
      for(auto const& seg : ktrk.fitTraj().pieces()){
        auto cdist = fabs(-1.0/seg->omega() - seg->d0());
        if( cdist < minCenterRho_){
          retval = false;
          break;
        }
      }
    }
    return retval;
  }

  void CentralHelixFit::sampleFit(KKTRK& ktrk) const {
    auto const& ftraj = ktrk.fitTraj();
    double tbeg = ftraj.range().begin();
    for(auto const& surf : surfacess_to_sample_){
      // search for intersections with each surface from the begining
      double tstart = tbeg - sampletbuff_;
      bool goodinter(true);
      size_t max_inter = 100;
      size_t cur_inter = 0;

      // loop to find multiple intersections
      while(goodinter && cur_inter < max_inter) {
        cur_inter += 1;
        TimeRange irange(tstart,std::max(ftraj.range().end(),tstart)+sampletbuff_);
        auto surfinter = KinKal::intersect(ftraj,*surf.second,irange,intertol_);
        goodinter = surfinter.onsurface_ && ( (! sampleinbounds_) || surfinter.inbounds_ ) && ( (!sampleinrange_) || irange.inRange(surfinter.time_));
        if(goodinter) {
          // save the intersection information
          ktrk.addIntersection(surf.first,surfinter);
          // update for the next intersection
          double epsilon = intertol_/ftraj.speed(surfinter.time_);
          tstart = surfinter.time_ + epsilon;// move psst existing intersection to avoid repeating
        }
      }
    }
  }

} // mu2e
using mu2e::CentralHelixFit;
DEFINE_ART_MODULE(CentralHelixFit)
