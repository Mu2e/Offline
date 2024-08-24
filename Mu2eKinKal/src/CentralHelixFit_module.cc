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
// C++
#include <iostream>
#include <fstream>
#include <string>
#include <functional>
#include <vector>
#include <memory>

using KTRAJ= KinKal::CentralHelix; // this must come before HelixFit
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
    fhicl::Atom<art::InputTag> timeClusterCollection         {Name("TimeClusterCollection"),     Comment("Seed fit collections to be processed ") };
    fhicl::Atom<art::InputTag> panelHitCollection         {Name("PanelHitCollection"),     Comment("Seed fit collections to be processed ") };
    fhicl::OptionalAtom<double> fixedBField { Name("ConstantBField"), Comment("Constant BField value") };
    fhicl::Atom<double> hitSep { Name("HitSeparation"), Comment("Seed with hits this far apart") };
    fhicl::Atom<double> seedMom { Name("SeedMomentum"), Comment("Seed momentum") };
    fhicl::Atom<int> seedCharge { Name("SeedCharge"), Comment("Seed charge") };
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
      TrkFitFlag fitflag_;
      // parameter-specific functions that need to be overridden in subclasses
      // data payload
      art::ProductToken<TimeClusterCollection> tccol_T_;
      art::ProductToken<ComboHitCollection> chcol_T_;
      art::ProductToken<ComboHitCollection> phcol_T_;
      art::ProductToken<CaloClusterCollection> cccol_T_;
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
      std::unique_ptr<KinKal::BFieldMap> kkbf_;
      Config config_; // initial fit configuration object
      Config exconfig_; // extension configuration object
      bool fixedfield_; //
      double hitSep_;
      double seedMom_;
      int seedCharge_;
    };

  CentralHelixFit::CentralHelixFit(const Parameters& settings) : art::EDProducer{settings},
    fitflag_(TrkFitFlag::KKCentralHelix),
    tccol_T_(consumes<TimeClusterCollection>(settings().modSettings().timeClusterCollection())),
    chcol_T_(consumes<ComboHitCollection>(settings().modSettings().comboHitCollection())),
    phcol_T_(consumes<ComboHitCollection>(settings().modSettings().panelHitCollection())),
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
    hitSep_(settings().modSettings().hitSep()),
    seedMom_(settings().modSettings().seedMom()),
    seedCharge_(settings().modSettings().seedCharge())
    {
      // collection handling
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
      if(print_ > 0) std::cout << "Fit " << config_ << "Extension " << exconfig_;
      double bz(0.0);
      if(settings().modSettings().fixedBField(bz)){
        fixedfield_ = true;
        kkbf_ = std::move(std::make_unique<KKConstantBField>(VEC3(0.0,0.0,bz)));
      }
    }

  void CentralHelixFit::beginRun(art::Run& run) {
    // setup things that rely on data related to beginRun
    auto const& ptable = GlobalConstantsHandle<ParticleDataList>();
    mass_ = ptable->particle(fpart_).mass();
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
    auto ph_H = event.getValidHandle<ComboHitCollection>(phcol_T_);
    auto cc_H = event.getValidHandle<CaloClusterCollection>(cccol_T_);
    auto const& chcol = *ch_H;
    auto const& phcol = *ph_H;
    // create output
    unique_ptr<KKTRKCOL> kktrkcol(new KKTRKCOL );
    unique_ptr<KalSeedCollection> kkseedcol(new KalSeedCollection );
//    unique_ptr<KalHelixAssns> kkseedassns(new KalHelixAssns());
    //auto KalSeedCollectionPID = event.getProductID<KalSeedCollection>();
    //auto KalSeedCollectionGetter = event.productGetter(KalSeedCollectionPID);

    auto  const& tcH = event.getValidHandle<TimeClusterCollection>(tccol_T_);
    const TimeClusterCollection& tccol(*tcH);

    for (size_t index=0;index<tccol.size();++index){
      const auto& tclust = tccol[index];
      ComboHitCollection combohits;
      combohits.setParent(phcol.parent());
      for (size_t i=0;i<tclust.hits().size();i++){
        combohits.push_back(phcol[tclust.hits()[i]]);
      }

      auto zcent = Mu2eKinKal::zMid(combohits);
      // take the magnetic field at the helix center as nominal
      auto bnom = kkbf_->fieldVect(VEC3(0,0,zcent));

      double bz = bnom.Z();
      auto trange = Mu2eKinKal::timeBounds(combohits);
      auto fitpart = fpart_;


      XYZVectorF best_seedpos, best_seeddir;
      double best_chi2 = 9e9;
      for (size_t i=0;i<combohits.size();i++){
        auto seedpos = combohits[i].pos();
        for (size_t j=i+1;j<combohits.size();j++){
          if ((combohits[j].pos()-combohits[i].pos()).Mag2() < hitSep_)
            continue;
          auto seeddir = (combohits[j].pos()-combohits[i].pos()).Unit();
          double chi2 = 0;
          for (size_t k=0;k<combohits.size();k++){
            if (k == i || k == j)
              continue;
            auto tohit = combohits[k].pos()-seedpos;
            auto poca = seedpos + seeddir*tohit.Dot(seeddir);
            auto delta = combohits[k].pos()-poca;
            double delta2 = delta.Mag2();
            delta = delta.Unit();

            double error2 = pow(delta.Dot(combohits[j].uDir())*combohits[j].uRes(),2)+pow(delta.Dot(combohits[j].vDir())*combohits[j].vRes(),2)+pow(delta.Dot(combohits[j].wDir())*combohits[j].wRes(),2);
            chi2 += delta2/error2;
          }
          if (chi2 < best_chi2){
            best_chi2 = chi2;
            best_seedpos = seedpos;
            best_seeddir = seeddir;
          }
        }
      }

      if (best_seeddir.Y() > 0)
        best_seeddir *= -1;

      double t0 = tclust._t0.t0();
      auto trackpos = best_seedpos;
      auto trackmom = best_seeddir * seedMom_;

      auto temptraj = KTRAJ(KinKal::VEC4(trackpos.x(),trackpos.y(),trackpos.z(),t0),KinKal::MOM4(trackmom.x(),trackmom.y(),trackmom.z(),mass_), seedCharge_, bz);
      KinKal::Parameters kkpars(temptraj.params().parameters(),seedcov_);
      auto seedtraj = KTRAJ(kkpars,mass_,seedCharge_,bz,trange);

      // wrap the seed traj in a Piecewise traj: needed to satisfy PTOCA interface
      PTRAJ pseedtraj(seedtraj);

      // first, we need to unwind the combohits.  We use this also to find the time range
      StrawHitIndexCollection strawHitIdxs;
      auto chcolptr = combohits.fillStrawHitIndices(strawHitIdxs, StrawIdMask::uniquestraw);
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
        auto kktrk = make_unique<KKTRK>(config_,*kkbf_,seedtraj,fitpart,kkfit_.strawHitClusterer(),strawhits,strawxings,calohits);
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
          kkseedcol->push_back(kkfit_.createSeed(*kktrk,fitflag,*calo_h));
          // save (unpersistable) KKTrk in the event
          kktrkcol->push_back(kktrk.release());
        }

      } catch (std::exception const& error) {
      }
    }
    // put the output products into the event
    if(print_ > 0) std::cout << "Fitted " << kktrkcol->size() << " tracks from " << tccol.size() << " Seeds" << std::endl;
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
    // test that the spatial parameter covariances and values aren't crazy TODO
    return retval;
  }



} // mu2e
using mu2e::CentralHelixFit;
DEFINE_ART_MODULE(CentralHelixFit)
