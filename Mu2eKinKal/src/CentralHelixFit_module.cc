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
#include "Offline/KinKalGeom/inc/SurfaceMap.hh"
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
#include "Offline/Mu2eKinKal/inc/KKMaterial.hh"
#include "Offline/Mu2eKinKal/inc/KKStrawHit.hh"
#include "Offline/Mu2eKinKal/inc/KKStrawHitCluster.hh"
#include "Offline/Mu2eKinKal/inc/KKStrawXing.hh"
#include "Offline/Mu2eKinKal/inc/KKCaloHit.hh"
#include "Offline/Mu2eKinKal/inc/KKBField.hh"
#include "Offline/Mu2eKinKal/inc/KKConstantBField.hh"
#include "Offline/Mu2eKinKal/inc/KKFitUtilities.hh"
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
  using KKFIT = KKFit<KTRAJ>;
  using KinKal::VEC3;
  using KinKal::DMAT;
  using HPtr = art::Ptr<HelixSeed>;
  using CCPtr = art::Ptr<CaloCluster>;
  using CCHandle = art::ValidHandle<CaloClusterCollection>;
  using StrawHitIndexCollection = std::vector<StrawHitIndex>;

  using KKConfig = Mu2eKinKal::KinKalConfig;
  using KKFitConfig = Mu2eKinKal::KKFitConfig;
  using KKModuleConfig = Mu2eKinKal::KKModuleConfig;

  using KKMaterialConfig = KKMaterial::Config;
  using Name    = fhicl::Name;
  using Comment = fhicl::Comment;
  using KinKal::DVEC;

  // extend the generic module configuration as needed
  struct KKCHModuleConfig : KKModuleConfig {
    fhicl::Sequence<art::InputTag> seedCollections         {Name("CosmicTrackSeedCollections"),     Comment("Seed fit collections to be processed ") };
    fhicl::OptionalAtom<double> fixedBField { Name("ConstantBField"), Comment("Constant BField value") };
    fhicl::Atom<double> seedMom { Name("SeedMomentum"), Comment("Seed momentum") };
    fhicl::Atom<int> seedCharge { Name("SeedCharge"), Comment("Seed charge MAGNITUDE, in electron charge units") };
    fhicl::Sequence<float> parconst { Name("ParameterConstraints"), Comment("External constraint on parameters to seed values (rms, various units)") };
    fhicl::Sequence<std::string> sampleSurfaces { Name("SampleSurfaces"), Comment("When creating the KalSeed, sample the fit at these surfaces") };
    fhicl::Atom<bool> sampleInRange { Name("SampleInRange"), Comment("Require sample times to be inside the fit trajectory time range") };
    fhicl::Atom<bool> sampleInBounds { Name("SampleInBounds"), Comment("Require sample intersection point be inside surface bounds (within tolerance)") };
    fhicl::Atom<float> sampleTol { Name("SampleTolerance"), Comment("Tolerance for sample surface intersections (mm)") };
    fhicl::Atom<float> sampleTBuff { Name("SampleTimeBuffer"), Comment("Time buffer for sample intersections (nsec)") };
    fhicl::Atom<bool> useFitCharge { Name("UseFitCharge"), Comment("Set the PDG particle according to the fit charge; otherwise reject fits that don't agree with the PDG particle charge") };
    fhicl::Atom<float> minCenterRho { Name("MinCenterRho"), Comment("Minimum transverse distance from the helix axis to the Z axis to consider the fit non-degenerate (mm)") };
};

  struct GlobalConfig {
    fhicl::Table<KKCHModuleConfig> modSettings { Name("ModuleSettings") };
    fhicl::Table<KKFitConfig> kkfitSettings { Name("KKFitSettings") };
    fhicl::Table<KKConfig> fitSettings { Name("FitSettings") };
    fhicl::Table<KKConfig> extSettings { Name("ExtensionSettings") };
    fhicl::Table<KKMaterialConfig> matSettings { Name("MaterialSettings") };
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
      void sampleFit(KKTRK const& kktrk,KalIntersectionCollection& inters) const;
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
      KKMaterial kkmat_; // material helper
      DMAT seedcov_; // seed covariance matrix
      double mass_; // particle mass
      int PDGcharge_; // PDG particle charge
      std::unique_ptr<KinKal::BFieldMap> kkbf_;
      Config config_; // initial fit configuration object
      Config exconfig_; // extension configuration object
      bool fixedfield_; //
      double seedMom_;
      int seedCharge_;
      double sampletol_; // surface intersection tolerance (mm)
      double sampletbuff_; // simple time buffer; replace this with extrapolation TODO
      bool useFitCharge_; // Set the PDG particle to agree with the fit charge
      double minCenterRho_; // min center distance to z axis
      bool sampleinrange_, sampleinbounds_; // require samples to be in range or on surface
      SurfaceMap::SurfacePairCollection sample_; // surfaces to sample the fit
      std::array<double,KinKal::NParams()> paramconstraints_;
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
    kkmat_(settings().matSettings()),
    config_(Mu2eKinKal::makeConfig(settings().fitSettings())),
    exconfig_(Mu2eKinKal::makeConfig(settings().extSettings())),
    fixedfield_(false),
    seedMom_(settings().modSettings().seedMom()),
    seedCharge_(settings().modSettings().seedCharge()),
    sampletol_(settings().modSettings().sampleTol()),
    sampletbuff_(settings().modSettings().sampleTBuff()),
    useFitCharge_(settings().modSettings().useFitCharge()),
    minCenterRho_(settings().modSettings().minCenterRho()),
    sampleinrange_(settings().modSettings().sampleInRange()),
    sampleinbounds_(settings().modSettings().sampleInBounds())
    {
      // collection handling
      for(const auto& cseedtag : settings().modSettings().seedCollections()) { cseedCols_.emplace_back(consumes<CosmicTrackSeedCollection>(cseedtag)); }
      produces<KKTRKCOL>();
      produces<KalSeedCollection>();
      //      produces<KalHelixAssns>();
      // build the initial seed covariance
      auto const& seederrors = settings().modSettings().seederrors();
      if(seederrors.size() != KinKal::NParams())
        throw cet::exception("RECO")<<"mu2e::CentralHelixFit:Seed error configuration error"<< endl;
      for(size_t ipar=0;ipar < seederrors.size(); ++ipar){
        seedcov_[ipar][ipar] = seederrors[ipar]*seederrors[ipar];
      }
      auto const& parerrors = settings().modSettings().parconst();
      for (size_t ipar=0;ipar<KinKal::NParams();ipar++)
        paramconstraints_[ipar] = parerrors[ipar];

      if(print_ > 0) std::cout << "Fit " << config_ << "Extension " << exconfig_;
      double bz(0.0);
      if(settings().modSettings().fixedBField(bz)){
        fixedfield_ = true;
        kkbf_ = std::move(std::make_unique<KKConstantBField>(VEC3(0.0,0.0,bz)));
      }
      SurfaceIdCollection ssids;
      for(auto const& sidname : settings().modSettings().sampleSurfaces()) {
        ssids.push_back(SurfaceId(sidname,-1)); // match all elements
      }
      // translate the sample and extend surface names to actual surfaces using the SurfaceMap.  This should come from the
      // geometry service eventually, TODO
      SurfaceMap smap;
      smap.surfaces(ssids,sample_);
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
  }

  void CentralHelixFit::produce(art::Event& event ) {
    GeomHandle<Calorimeter> calo_h;
    // find current proditions
    auto const& strawresponse = strawResponse_h_.getPtr(event.id());
    auto const& tracker = alignedTracker_h_.getPtr(event.id()).get();
    // find input hits
    auto ch_H = event.getValidHandle<ComboHitCollection>(chcol_T_);
    auto cc_H = event.getValidHandle<CaloClusterCollection>(cccol_T_);
    auto const& chcol = *ch_H;
    // create output
    unique_ptr<KKTRKCOL> kktrkcol(new KKTRKCOL );
    unique_ptr<KalSeedCollection> kkseedcol(new KalSeedCollection );

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
        kkfit_.makeStrawHits(*tracker, *strawresponse, *kkbf_, kkmat_.strawMaterial(), pseedtraj, chcol, strawHitIdxs, strawhits, strawxings);
        // optionally (and if present) add the CaloCluster as a constraint
        // verify the cluster looks physically reasonable before adding it TODO!  Or, let the KKCaloHit updater do it TODO
        KKCALOHITCOL calohits;
        //FIXME    if (kkfit_.useCalo() && hseed.caloCluster().isNonnull())kkfit_.makeCaloHit(hseed.caloCluster(),*calo_h, pseedtraj, calohits);
        // extend the seed range given the hits and xings

        try {
          seedtraj.range() = kkfit_.range(strawhits,calohits,strawxings);

          // create and fit the track
          auto kktrk = make_unique<KKTRK>(config_,*kkbf_,seedtraj,fitpart,kkfit_.strawHitClusterer(),strawhits,strawxings,calohits,paramconstraints_);
          // Check the fit
          auto goodfit = goodFit(*kktrk);
          // if we have an extension schedule, extend.
          if(goodfit && exconfig_.schedule().size() > 0) {
            //  std::cout << "EXTENDING TRACK " << event.id() << " " << index << std::endl;
            kkfit_.extendTrack(exconfig_,*kkbf_, *tracker,*strawresponse, kkmat_.strawMaterial(), chcol, *calo_h, cc_H, *kktrk );
            goodfit = goodFit(*kktrk);
          }

          if(print_>1)kktrk->printFit(std::cout,print_);
          if(goodfit || saveall_){
            TrkFitFlag fitflag;//(hptr->status());
            fitflag.merge(fitflag_);
            if(goodfit)
              fitflag.merge(TrkFitFlag::FitOK);
            else
              fitflag.clear(TrkFitFlag::FitOK);
            // compare charge after the fit; either adjust or skip
            auto const& t0seg = kktrk->fitTraj().nearestPiece(kktrk->fitTraj().t0());
            // check that the t0 segment is non-degenerate; skip tracks that are. This must be done in helix-local coordinates
            auto g2l = Rotation3D(AxisAngle(VEC3(sin(t0seg.bnom().Phi()),-cos(t0seg.bnom().Phi()),0.0),t0seg.bnom().Theta()));
            auto cpos = g2l(t0seg.center(t0seg.range().mid()));
            auto cdist = cpos.Rho();
            if( cdist > minCenterRho_){
              double t0charge = t0seg.charge();
              if(t0charge*PDGcharge_> 0 || useFitCharge_){
                // flip the PDG particle assignment charge if required
                if(t0charge*PDGcharge_ < 0)kktrk->reverseCharge();
                auto kkseed = kkfit_.createSeed(*kktrk,fitflag,*calo_h);
                sampleFit(*kktrk,kkseed._inters);
                kkseedcol->push_back(kkseed);
                // save (unpersistable) KKTrk in the event
                kktrkcol->push_back(kktrk.release());
              }
            }
          }

        } catch (std::invalid_argument const& error) {
          if(print_ > 0) std::cout << "CentralHelixFit Error " << error.what() << std::endl;
        }
      }
    }
    // put the output products into the event
    if(print_ > 0) std::cout << "Fitted " << kktrkcol->size() << " tracks from " << nseed << " Seeds" << std::endl;
    event.put(move(kktrkcol));
    event.put(move(kkseedcol));
  }

  bool CentralHelixFit::goodFit(KKTRK const& ktrk) const {
    // require physical consistency: fit can succeed but the result can have changed charge or helicity
    bool retval = ktrk.fitStatus().usable();
    // also check that the fit is inside the physical detector volume.  Test where the StrawHits are
    if(retval){
      for(auto const& shptr : ktrk.strawHits()) {
        if(shptr->active() && !Mu2eKinKal::inDetector(ktrk.fitTraj().position3(shptr->time()))){
          retval = false;
          break;
        }
      }
    }
    return retval;
  }

  void CentralHelixFit::sampleFit(KKTRK const& kktrk,KalIntersectionCollection& inters) const {
    auto const& ftraj = kktrk.fitTraj();
    double tbeg = ftraj.range().begin();
    static const double epsilon(1.0e-3);
    for(auto const& surf : sample_){
      // search for intersections with each surface from the begining
      double tstart = tbeg - sampletbuff_;
      bool hasinter(true);
      size_t max_iter = 1000;
      size_t cur_iter = 0;

      // loop to find multiple intersections
      while(hasinter) {
        if (cur_iter > max_iter)
          break;
        cur_iter += 1;

        TimeRange irange(tstart,std::max(ftraj.range().end(),tstart)+sampletbuff_);
        auto surfinter = KinKal::intersect(ftraj,*surf.second,irange,sampletol_);
        hasinter = surfinter.onsurface_ && ( (! sampleinbounds_) || surfinter.inbounds_ ) && ( (!sampleinrange_) || irange.inRange(surfinter.time_));
        if(hasinter) {
          // save the intersection information
          auto const& ktraj = ftraj.nearestPiece(surfinter.time_);
          inters.emplace_back(ktraj.stateEstimate(surfinter.time_),XYZVectorF(ktraj.bnom()),surf.first,surfinter);
          // update for the next intersection
          tstart = surfinter.time_ + epsilon;// move psst existing intersection to avoid repeating
        }
      }
    }
    // record other intersections saved in the track
    for(auto const& interpair : kktrk.intersections()) {
      auto const& sid = std::get<0>(interpair);
      auto const& inter = std::get<1>(interpair);
      auto const& ktraj = ftraj.nearestPiece(inter.time_);
      inters.emplace_back(ktraj.stateEstimate(inter.time_),XYZVectorF(ktraj.bnom()),sid,inter);
    }
    // sort by time TODO
  }

} // mu2e
using mu2e::CentralHelixFit;
DEFINE_ART_MODULE(CentralHelixFit)
