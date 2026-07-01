//
// KinKal fit module using the KinematicLine parameterset
//
// Original author S. Middleton (Caltech) 2021
//
// framework
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/OptionalTable.h"
#include "fhiclcpp/types/Tuple.h"
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
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/Mu2eUtilities/inc/CosmicTrackUtils.hh"
// data
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/DataProducts/inc/Helicity.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/CosmicTrackSeed.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/KKLine.hh"
#include "Offline/DataProducts/inc/SurfaceId.hh"
#include "Offline/KinKalGeom/inc/KinKalGeom.hh"
// MC (TruthSeedDiag only)
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
// KinKal
#include "KinKal/Fit/Track.hh"
#include "KinKal/Fit/Config.hh"
#include "KinKal/Trajectory/KinematicLine.hh"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "KinKal/Trajectory/PiecewiseClosestApproach.hh"
#include "KinKal/General/Parameters.hh"
// Mu2eKinKal
#include "Offline/Mu2eKinKal/inc/KKFit.hh"
#include "Offline/Mu2eKinKal/inc/KKFitSettings.hh"
#include "Offline/KinKalGeom/inc/KKMaterial.hh"
#include "Offline/Mu2eKinKal/inc/KKStrawHit.hh"
#include "Offline/Mu2eKinKal/inc/KKBField.hh"
#include "Offline/Mu2eKinKal/inc/KKFitUtilities.hh"
#include "Offline/Mu2eKinKal/inc/KKShellXing.hh"
#include "Offline/Mu2eKinKal/inc/KKExtrap.hh"
// root
#include "TH1F.h"
#include "TTree.h"
// C++
#include <iostream>
#include <fstream>
#include <string>
#include <functional>
#include <vector>
#include <memory>
#include <limits>
using namespace std;
//using namespace KinKal;
namespace mu2e {
  using KTRAJ= KinKal::KinematicLine;
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
  using MEAS = KinKal::Hit<KTRAJ>;
  using MEASPTR = std::shared_ptr<MEAS>;
  using MEASCOL = std::vector<MEASPTR>;
  using KKFIT = KKFit<KTRAJ>;
  using EXING = KinKal::ElementXing<KTRAJ>;
  using EXINGPTR = std::shared_ptr<EXING>;
  using EXINGCOL = std::vector<EXINGPTR>;
  using KKFIT = mu2e::KKFit<KTRAJ>;
  using KinKal::DVEC;
  using KinKal::VEC3;
  using KinKal::TimeRange;
  using KinKal::DMAT;
  using KinKal::Status;
  using HPtr = art::Ptr<CosmicTrackSeed>;
  using CCPtr = art::Ptr<CaloCluster>;
  using CCHandle = art::Handle<CaloClusterCollection>;
  using StrawHitIndexCollection = std::vector<StrawHitIndex>;
  using KKCRVXING = KKShellXing<KTRAJ,KinKal::Rectangle>;
  using KKCRVXINGPTR = std::shared_ptr<KKCRVXING>;
  using KKCRVXINGCOL = std::vector<KKCRVXINGPTR>;

  using KKConfig = Mu2eKinKal::KinKalConfig;
  using KKFitConfig = Mu2eKinKal::KKFitConfig;
  using KKModuleConfig = Mu2eKinKal::KKModuleConfig;

  class KinematicLineFit : public art::EDProducer {
    using Name    = fhicl::Name;
    using Comment = fhicl::Comment;
  // extend the generic module configuration as needed
    struct KKLineModuleConfig : public KKModuleConfig {
      fhicl::Sequence<art::InputTag> seedCollections         {Name("CosmicTrackSeedCollections"),     Comment("Seed fit collections to be processed ") };
      fhicl::Atom<float> seedmom { Name("SeedMomentum"), Comment("Initial momentum value")};
      fhicl::Sequence<double> paramconstraints { Name("ParameterConstraints"), Comment("Sigma of direct gaussian constraints on each parameter (0=no constraint)")};
      fhicl::Sequence<std::string> sampleSurfaces { Name("SampleSurfaces"), Comment("When creating the KalSeed, sample the fit at these surfaces") };
      fhicl::Atom<bool> sampleInRange { Name("SampleInRange"), Comment("Require sample times to be inside the fit trajectory time range") };
      fhicl::Atom<bool> sampleInBounds { Name("SampleInBounds"), Comment("Require sample intersection point be inside surface bounds (within tolerance)") };
      fhicl::Atom<float> interTol { Name("IntersectionTolerance"), Comment("Tolerance for surface intersections (mm)") };
      fhicl::Atom<float> sampleTBuff { Name("SampleTimeBuffer"), Comment("Time buffer for sample intersections (nsec)") };
      fhicl::Atom<bool> truthSeedDiag { Name("TruthSeedDiag"), Comment("Diagnostic: also extrapolate a copy of each fit re-seeded with the true muon state, isolating extrapolation error from fit error. Writes KalSeeds under instance 'TruthSeed'"), false };
      fhicl::OptionalAtom<art::InputTag> vdStepCollection { Name("VDStepPointMCs"), Comment("Virtual-detector StepPointMCs providing the true muon state at the tracker (required iff TruthSeedDiag)") };
      fhicl::Sequence<double> truthSeedParamConstraints {
        Name("TruthSeedParameterConstraints"),
        Comment("TruthSeedDiag ParameterHit RMS constraints for KinematicLine parameters d0, phi0, z0, cost, t0, mom"),
        std::vector<double>{1.0e-3, 1.0e-6, 1.0e-3, 1.0e-6, 1.0e-3, 1.0e-3}
      };
    };

    struct GlobalConfig {
      fhicl::Table<KKLineModuleConfig> modSettings { Name("ModuleSettings") };
      fhicl::Table<KKFitConfig> mu2eSettings { Name("KKFitSettings") };
      fhicl::Table<KKConfig> fitSettings { Name("FitSettings") };
      fhicl::Table<KKConfig> extSettings { Name("ExtensionSettings") };
      fhicl::OptionalTable<KKExtrapConfig> extrapSettings { Name("ExtrapolationSettings") };
      fhicl::OptionalAtom<std::string> fitDirection { Name("FitDirection"), Comment("Particle direction to fit, either \"upstream\" or \"downstream\"")};
    };

    public:
    using Parameters = art::EDProducer::Table<GlobalConfig>;
    explicit KinematicLineFit(const Parameters& settings);
    virtual ~KinematicLineFit();
    void beginRun(art::Run& run) override;
    void produce(art::Event& event) override;
    private:
    // utility functions
    KTRAJ makeSeedTraj(CosmicTrackSeed const& hseed) const;
    bool goodFit(KKTRK const& ktrk) const;
    void sampleFit(KKTRK& kktrk) const;
    void extrapolate(KKTRK& ktrk) const;
    // data payload
    std::vector<art::ProductToken<CosmicTrackSeedCollection>> seedCols_;
    art::ProductToken<ComboHitCollection> chcol_T_;
    art::ProductToken<CaloClusterCollection> cccol_T_;
    TrkFitFlag goodline_;
    bool saveall_;
    ProditionsHandle<StrawResponse> strawResponse_h_;
    ProditionsHandle<Tracker> alignedTracker_h_;
    int print_;
    float seedmom_;
    std::vector<double> paramconstraints_;
    PDGCode::type fpart_;
    TrkFitDirection fdir_;
    KKFIT kkfit_; // fit helper
    DMAT seedcov_; // seed covariance matrix
    bool constraining_;
    double mass_; // particle mass
    int charge_; // particle charge
    std::unique_ptr<KKBField> kkbf_;
    double intertol_; // surface intersection tolerance (mm)
    double sampletbuff_; // simple time buffer; replace this with extrapolation TODO
    bool sampleinrange_, sampleinbounds_; // require samples to be in range or on surface
    std::unique_ptr<KKExtrap> extrap_; // extrapolation helper
    SurfaceIdCollection ssids_;
    KinKalGeom::SurfacePairCollection surfacess_to_sample_; // surfaces to sample the fit
    Config config_; // initial fit configuration object
    Config exconfig_; // extension configuration object
    bool truthSeedDiag_; // diagnostic: extrapolate a truth-seeded copy of each fit
    art::ProductToken<StepPointMCCollection> vdcol_T_; // VD StepPointMCs (truth seed)
    std::vector<double> truthSeedParamConstraints_; // ParameterHit RMS constraints for truth seeding
  };

  KinematicLineFit::KinematicLineFit(const Parameters& settings) : art::EDProducer{settings},
    chcol_T_(consumes<ComboHitCollection>(settings().modSettings().comboHitCollection())),
    cccol_T_(mayConsume<CaloClusterCollection>(settings().modSettings().caloClusterCollection())),
    goodline_(settings().modSettings().seedFlags()),
    saveall_(settings().modSettings().saveAll()),
    print_(settings().modSettings().printLevel()),
    seedmom_(settings().modSettings().seedmom()),
    paramconstraints_(settings().modSettings().paramconstraints()),
    fpart_(static_cast<PDGCode::type>(settings().modSettings().fitParticle())),
    kkfit_(settings().mu2eSettings()),
    intertol_(settings().modSettings().interTol()),
    sampletbuff_(settings().modSettings().sampleTBuff()),
    sampleinrange_(settings().modSettings().sampleInRange()),
    sampleinbounds_(settings().modSettings().sampleInBounds()),
    config_(Mu2eKinKal::makeConfig(settings().fitSettings())),
    exconfig_(Mu2eKinKal::makeConfig(settings().extSettings())),
    truthSeedDiag_(settings().modSettings().truthSeedDiag()),
    truthSeedParamConstraints_(settings().modSettings().truthSeedParamConstraints())
    {
      std::string fdir;
      if(settings().fitDirection(fdir))fdir_ = fdir;
      // collection handling
      for(const auto& seedtag : settings().modSettings().seedCollections()) { seedCols_.emplace_back(consumes<CosmicTrackSeedCollection>(seedtag)); }
      produces<KKTRKCOL>();
      produces<KalSeedCollection>();
      produces<KalLineAssns>();
      // TruthSeedDiag: a parallel truth-seeded extrapolation, written under instance "TruthSeed"
      if(truthSeedDiag_){
        art::InputTag vdtag;
        if(!settings().modSettings().vdStepCollection(vdtag))
          throw cet::exception("RECO") << "mu2e::KinematicLineFit: TruthSeedDiag requires VDStepPointMCs" << endl;
        if(truthSeedParamConstraints_.size() != KinKal::NParams())
          throw cet::exception("RECO") << "mu2e::KinematicLineFit: TruthSeedParameterConstraints must have "
            << KinKal::NParams() << " entries" << endl;
        for(auto const& sigma : truthSeedParamConstraints_){
          if(sigma <= 0.0)
            throw cet::exception("RECO") << "mu2e::KinematicLineFit: TruthSeedParameterConstraints entries must be positive" << endl;
        }
        vdcol_T_ = consumes<StepPointMCCollection>(vdtag);
        produces<KKTRKCOL>("TruthSeed");
        produces<KalSeedCollection>("TruthSeed");
      }
      // build the initial seed covariance
      auto const& seederrors = settings().modSettings().seederrors();
      if(seederrors.size() != KinKal::NParams())
        throw cet::exception("RECO")<<"mu2e::KinematicLineFit:Seed error configuration error"<< endl;
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
        throw cet::exception("RECO")<<"mu2e::KinematicLineFit: Parameter constraint configuration error"<< endl;
      }
      for(auto const& sidname : settings().modSettings().sampleSurfaces()) {
        ssids_.push_back(SurfaceId(sidname,-1)); // match all elements
      }
      // setup extrapolation
      if(settings().extrapSettings())extrap_ = make_unique<KKExtrap>(*settings().extrapSettings());
      if(print_ > 0) std::cout << config_;
    }

  KinematicLineFit::~KinematicLineFit(){}

  void KinematicLineFit::beginRun(art::Run& run) {
    // setup particle parameters
    auto const& ptable = GlobalConstantsHandle<ParticleDataList>();
    mass_ = ptable->particle(fpart_).mass();
    charge_ = static_cast<int>(ptable->particle(fpart_).charge());
    // create KKBField
    GeomHandle<BFieldManager> bfmgr;
    GeomHandle<DetectorSystem> det;
    kkbf_ = std::make_unique<KKBField>(*bfmgr,*det);
    // translate the sample surface names to actual surfaces using the KinKalGeom. This must be done after construction as the KKGeom object now comes from GeometryService
    GeomHandle<mu2e::KinKalGeom> kkg_h;
    kkg_h->surfaces(ssids_,surfacess_to_sample_);
  }

  void KinematicLineFit::produce(art::Event& event ) {
    GeomHandle<mu2e::Calorimeter> calo_h;
    GeomHandle<mu2e::Tracker> nominalTracker_h;

    // find current proditions
    auto const& strawresponse = strawResponse_h_.getPtr(event.id());
    auto const& tracker = alignedTracker_h_.getPtr(event.id()).get();
    // find input hits
    auto ch_H = event.getValidHandle<ComboHitCollection>(chcol_T_);
    auto cc_H = event.getHandle<CaloClusterCollection>(cccol_T_);
    auto const& chcol = *ch_H;
    // create output
    unique_ptr<KKTRKCOL> kktrkcol(new KKTRKCOL );
    unique_ptr<KalSeedCollection> kkseedcol(new KalSeedCollection ); //Needs to return a KalSeed
    unique_ptr<KalLineAssns> kkseedassns(new KalLineAssns());
    // TruthSeedDiag parallel outputs + truth inputs
    unique_ptr<KKTRKCOL> kktrkcol_truth(new KKTRKCOL );
    unique_ptr<KalSeedCollection> kkseedcol_truth(new KalSeedCollection );
    StepPointMCCollection const* vdsteps = nullptr;
    GeomHandle<DetectorSystem> det;
    if(truthSeedDiag_) vdsteps = event.getValidHandle<StepPointMCCollection>(vdcol_T_).product();
    auto KalSeedCollectionPID = event.getProductID<KalSeedCollection>();
    auto KalSeedCollectionGetter = event.productGetter(KalSeedCollectionPID);
    // find the track seed collections
    unsigned nseed(0);
    for (auto const& hseedtag : seedCols_) {
      auto const& hseedcol_h = event.getValidHandle<CosmicTrackSeedCollection>(hseedtag);
      auto const& hseedcol = *hseedcol_h;
      nseed += hseedcol.size();
      // loop over the seeds
      for(size_t iseed=0; iseed < hseedcol.size(); ++iseed) {
        auto const& hseed = hseedcol[iseed];
        auto hptr = HPtr(hseedcol_h,iseed);
        // check helicity.  The test on the charge and helicity
        if(hseed.status().hasAllProperties(goodline_) ){
          // construt the seed trajectory
          KTRAJ seedtraj = makeSeedTraj(hseed);
          // wrap the seed traj in a Piecewise traj: needed to satisfy PTOCA interface
          PTRAJ pseedtraj(seedtraj);
          // first, we need to unwind the combohits.  We use this also to find the time range
          StrawHitIndexCollection strawHitIdxs;
          auto chcolptr = hseed.hits().fillStrawHitIndices(strawHitIdxs, StrawIdMask::uniquestraw);
//          if(chcolptr != &chcol)
//            throw cet::exception("RECO")<<"mu2e::KinematicLineFit: inconsistent ComboHitCollection" << std::endl;
          // next, build straw hits and materials from these
          KKSTRAWHITCOL strawhits;
          KKSTRAWXINGCOL strawxings;
          strawhits.reserve(strawHitIdxs.size());
          strawxings.reserve(strawHitIdxs.size());
          kkfit_.makeStrawHits(*tracker, *strawresponse, *kkbf_, pseedtraj, *chcolptr, strawHitIdxs, strawhits, strawxings);

          //here
          KKCALOHITCOL calohits;
          //if (kkfit_.useCalo()) kkfit_.makeCaloHit(hptr->caloCluster(),*calo_h, pseedtraj, calohits); --> CosmicTrackSeed has no CaloClusters....

          if(print_ > 2){
            for(auto const& strawhit : strawhits) strawhit->print(std::cout,2);
            for(auto const& calohit : calohits) calohit->print(std::cout,2);
            for(auto const& strawxing :strawxings) strawxing->print(std::cout,2);
          }
          // set the seed range given the hit TPOCA values
          seedtraj.range() = kkfit_.range(strawhits,calohits, strawxings);

          // create parameter constraint, use seedtrajectory parameters
          PARAMHITCOL paramhits;
          if (constraining_){
            kkfit_.makeSeedParamHit(seedtraj,paramconstraints_,paramhits);
          }

          if(print_ > 0){
            //std::cout << "Seed line parameters " << hseed.track() << std::endl;
            seedtraj.print(std::cout,print_);
          }
          // create and fit the track
          auto kktrk = make_unique<KKTRK>(config_,*kkbf_,seedtraj,fpart_,kkfit_.strawHitClusterer(),strawhits,strawxings,calohits,paramhits);
          auto goodfit = goodFit(*kktrk);
          if(goodfit && exconfig_.schedule().size() > 0){
            kkfit_.extendTrack(exconfig_,*kkbf_, *tracker,*strawresponse, chcol, *calo_h, cc_H, *kktrk );
          }
          goodfit = goodFit(*kktrk);
          // [TruthSeedDiag] Build a parallel track seeded with the TRUE muon state at the tracker
          // (nearest VD crossing), constrained to it by ParameterHits, and extrapolate that -- via the
          // PUBLIC KinKal reconstruction interface (no fitTrajMutable). The truth-seeded CRV residual
          // isolates the extrapolation (material/scattering) error from the fit's contribution.
          if(truthSeedDiag_ && goodfit && extrap_ && vdsteps != nullptr){
            double tt0 = kktrk->fitTraj().t0();
            auto refpos = kktrk->fitTraj().position3(tt0);
            double bestd2 = std::numeric_limits<double>::max();
            VEC3 tpos, tmom; double ttime = 0.0; int tcharge = 0; bool found = false;
            for(auto const& vds : *vdsteps){ // nearest true muon VD crossing to the fit t0 point
              auto const& sp = vds.simParticle();
              if(sp.isNull() || std::abs(static_cast<int>(sp->pdgId())) != PDGCode::mu_minus) continue;
              auto const hd = det->toDetector(vds.position());
              VEC3 dpos(hd.x(),hd.y(),hd.z());
              double d2 = (dpos-refpos).Mag2();
              if(d2 < bestd2){
                bestd2 = d2; found = true; ttime = vds.time();
                tpos = dpos; tmom = VEC3(vds.momentum().x(),vds.momentum().y(),vds.momentum().z());
                tcharge = (static_cast<int>(sp->pdgId()) > 0) ? -1 : 1; // mu-(+13)->-1, mu+(-13)->+1
              }
            }
            if(found){
              try {
                // field-off straight line: same non-zero nominal bfield (VEC3) as makeSeedTraj
                VEC3 tbnom(0.0,0.0,0.001);
                KTRAJ truthtraj(KinKal::VEC4(tpos.X(),tpos.Y(),tpos.Z(),ttime),
                                KinKal::MOM4(tmom.X(),tmom.Y(),tmom.Z(),mass_),
                                tcharge, tbnom, kktrk->fitTraj().range());
                truthtraj.params() = KinKal::Parameters(truthtraj.params().parameters(), seedcov_);
                PTRAJ truthpt(truthtraj);
                // constrain the line to the truth params via ParameterHits (public reco interface). A
                // single piece is a zero-length fit range (lowNDOF), so pin TWO hits a sliver apart.
                PARAMHITCOL truthParamHits;
                PARAMHIT::PMASK truthMask; truthMask.fill(true);
                auto addTruthHit = [&](double hitTime, double covarianceScale){
                  auto cparams = truthpt.nearestPiece(hitTime).params();
                  for(size_t ipar = 0; ipar < KinKal::NParams(); ++ipar){
                    for(size_t jpar = 0; jpar < KinKal::NParams(); ++jpar) cparams.covariance()[ipar][jpar] = 0.0;
                    auto const sigma = truthSeedParamConstraints_.at(ipar);
                    cparams.covariance()[ipar][ipar] = covarianceScale*sigma*sigma;
                  }
                  truthParamHits.push_back(std::make_shared<PARAMHIT>(hitTime, truthpt, cparams, truthMask));
                };
                constexpr double truthHitDt = 1.0e-3;
                // Pin the constraints at an IN-RANGE time. The KKLine fit is field-off, so
                // kktrk->domains() is empty and Track::fit's detrange collapses to just these hit
                // times (Track.hh:202-214 only extend detrange to surviving domains). The VD
                // crossing time `ttime` can fall outside the fit trajectory's time range; pinning
                // the hits there leaves the single truth piece (range = fitTraj range) not
                // overlapping detrange, so PiecewiseTrajectory::setRange(trim) erases every piece
                // and then dereferences the empty deque -> SIGSEGV. The single-piece KinematicLine
                // params are time-independent, so clamping the constraint time into the fit range
                // pins exactly the same truth params. (CentralHelix avoids this via field domains +
                // in-range piece-midpoint sampling.)
                auto const& frange = kktrk->fitTraj().range();
                double tctr = std::max(frange.begin() + truthHitDt, std::min(ttime, frange.end() - truthHitDt));
                addTruthHit(tctr - 0.5*truthHitDt, 2.0);
                addTruthHit(tctr + 0.5*truthHitDt, 2.0);
                KKTRK::DOMAINCOL truthDomains;
                for(auto const& domain : kktrk->domains()) truthDomains.emplace(std::make_shared<KinKal::Domain>(*domain));
                KKTRK::PKTRAJPTR truthTraj = std::make_unique<PTRAJ>(truthpt);
                KKSTRAWHITCOL nostrawhits; KKSTRAWXINGCOL nostrawxings; KKCALOHITCOL nocalohits;
                // one iteration, no convergence churn: the ParameterHits already pin the params
                auto truthConfig = config_;
                truthConfig.maxniter_ = 1;
                truthConfig.divdchisq_ = std::numeric_limits<double>::max();
                truthConfig.pdchisq_ = std::numeric_limits<double>::max();
                truthConfig.divgap_ = std::numeric_limits<double>::max();
                auto kktrk_truth = make_unique<KKTRK>(truthConfig,*kkbf_,fpart_,truthTraj,
                    nostrawhits,nostrawxings,nocalohits,truthParamHits,truthDomains);
                if(kktrk_truth->fitStatus().usable()){
                  extrap_->extrapolate(*kktrk_truth);
                  TrkFitFlag tflag(hptr->status()); tflag.merge(TrkFitFlag::KKLine);
                  auto truthKSeed = kkfit_.createSeed(*kktrk_truth,tflag,*calo_h,*nominalTracker_h);
                  // carry the reco fit's hit metadata onto the truth KalSeed (the truth-constrained
                  // fit has no detector hits, but EventNtuple consumers expect the hit list)
                  auto recoMetadata = kkfit_.createSeed(*kktrk,tflag,*calo_h,*nominalTracker_h);
                  truthKSeed._hits = std::move(recoMetadata._hits);
                  truthKSeed._hitcalibs = std::move(recoMetadata._hitcalibs);
                  truthKSeed._chit = std::move(recoMetadata._chit);
                  kkseedcol_truth->push_back(std::move(truthKSeed));
                  kktrkcol_truth->push_back(kktrk_truth.release());
                } else if(print_ > 0){
                  std::cout << "KinematicLineFit TruthSeed ParameterHit fit unusable: " << kktrk_truth->fitStatus() << std::endl;
                }
              } catch (std::exception const& error) {
                if(print_ > 0) std::cout << "KinematicLineFit TruthSeed Error " << error.what() << std::endl;
              }
            }
          }
          // extrapolate as required
          if(goodfit && extrap_) extrap_->extrapolate(*kktrk);
          bool save = goodFit(*kktrk);
          if(save || saveall_){
            TrkFitFlag fitflag(hptr->status());
            fitflag.merge(TrkFitFlag::KKLine);
            sampleFit(*kktrk);
            auto kkseed = kkfit_.createSeed(*kktrk,fitflag,*calo_h,*nominalTracker_h);
            kkseedcol->push_back(kkseed);
            kkseedcol->back()._status.merge(TrkFitFlag::KKLine);
            // fill assns with the cosmic seed
            auto kseedptr = art::Ptr<KalSeed>(KalSeedCollectionPID,kkseedcol->size()-1,KalSeedCollectionGetter);
            kkseedassns->addSingle(kseedptr,hptr);
            // save (unpersistable) KKTrk in the event
            kktrkcol->push_back(kktrk.release());
          }
        }
      }
    }
    // put the output products into the event
    if(print_ > 0) std::cout << "Fitted " << kktrkcol->size() << " tracks from " << nseed << " seeds" << std::endl;
    event.put(move(kktrkcol));
    event.put(move(kkseedcol));
    event.put(move(kkseedassns));
    if(truthSeedDiag_){
      event.put(move(kktrkcol_truth),"TruthSeed");
      event.put(move(kkseedcol_truth),"TruthSeed");
    }
  }

  KTRAJ KinematicLineFit::makeSeedTraj(CosmicTrackSeed const& hseed) const {
    //exctract CosmicTrack (contains parameters)
    VEC3 bnom(0.0,0.0,0.001);// non-zero value doesn't affect fit, but insures consistency with interfaces.
    KinKal::VEC4 pos(hseed._track.MinuitParams.A0, 0, hseed._track.MinuitParams.B0, hseed._t0._t0);
    XYZVectorF mom3(hseed._track.MinuitParams.A1, -1, hseed._track.MinuitParams.B1);
    mom3 = mom3.Unit()*seedmom_;
    if (fdir_ == TrkFitDirection::FitDirection::downstream){
      if (mom3.z() < 0)
        mom3 *= -1;
    }else if (fdir_ == TrkFitDirection::FitDirection::upstream){
      if (mom3.z() > 0)
        mom3 *= -1;
    }
    KinKal::MOM4 mom(mom3.x(),mom3.y(),mom3.z(),mass_);

    auto seedtraj = KTRAJ(pos,mom,charge_,bnom,Mu2eKinKal::timeBounds(hseed.hits()));
    seedtraj.params().covariance() = seedcov_;
    return seedtraj;
  }
  bool KinematicLineFit::goodFit(KKTRK const& ktrk) const {
    // require physical consistency: fit can succeed but the result can have changed charge or helicity
    return ktrk.fitStatus().usable();
  }

  void KinematicLineFit::sampleFit(KKTRK& kktrk) const {
    auto const& ftraj = kktrk.fitTraj();
    // need to extend range for now even if sampleinrange_ is false
    TimeRange extrange(ftraj.range().begin() - sampletbuff_,ftraj.range().end() + sampletbuff_);
    kktrk.extendTraj(extrange);
    double tbeg = ftraj.range().begin();

    for(auto const& surf : surfacess_to_sample_){
      // search for intersections with each surface from the begining
      double tstart = tbeg;
      bool goodinter(true);
      size_t max_inter = 100;
      size_t cur_inter = 0;

      // loop to find multiple intersections
      while(goodinter && cur_inter < max_inter) {
        cur_inter += 1;
        TimeRange irange(tstart,std::max(ftraj.range().end(),tstart));
        auto surfinter = KinKal::intersect(ftraj,*surf.second,irange,intertol_);
        goodinter = surfinter.onsurface_ && ( (! sampleinbounds_) || surfinter.inbounds_ ) && ( (!sampleinrange_) || irange.inRange(surfinter.time_));
        if(goodinter) {
          // save the intersection information
          kktrk.addIntersection(surf.first,surfinter);
          // move past existing intersection to avoid repeating
          double epsilon = intertol_/ftraj.speed(surfinter.time_);
          // update for the next intersection
          tstart = surfinter.time_ + epsilon;// move psst existing intersection to avoid repeating
        }
      }
    }
  }
}
DEFINE_ART_MODULE(mu2e::KinematicLineFit)
