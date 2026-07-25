//
//  Investigating relationship between fitted helices
//  and calorimeter clustes with a view to skipping photon hits
//  Iffat Zarif, Nikita Mazotov
//
// framework
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"


#include "fhiclcpp/types/Atom.h"

// mu2e
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"
#include "Offline/RecoDataProducts/inc/TriggerInfo.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/CaloCluster/inc/ClusterUtils.hh"

// MC
#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/PrimaryParticle.hh"
#include "Offline/MCDataProducts/inc/GenId.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/MCDataProducts/inc/CaloMCTruthAssns.hh"
#include "Offline/MCDataProducts/inc/CaloClusterMC.hh"


// ROOT
#include "TH1F.h"

// c++
#include <iostream>
#include <cmath>
#include <limits>

namespace mu2e {

  class PhotonFilter : public art::EDFilter {
  public:
    struct Config {
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<art::InputTag> CaloClusterCollection {
        Name("CaloClusterCollection"), Comment("InputTag for CaloClusterCollection")
      };

      fhicl::Atom<art::InputTag> HelixSeedCollection {
        Name("HelixSeedCollection"), Comment("InputTag for HelixSeedCollection")
      };

      fhicl::Atom<art::InputTag>  CaloClusterMC          {
        Name("CaloClusterMC"),           Comment("CaloCluster truth name")
      };
      fhicl::Atom<art::InputTag>  PrimaryParticle          {
        Name("PrimaryParticle"),           Comment("Primary Particle")
      };

      fhicl::Atom<int> DebugLevel {
        Name("DebugLevel"), Comment("Debug verbosity level (0=quiet, 1=verbose)"), 0
      };
      fhicl::Atom<double> MinDeltaR {
        Name("MinDeltaR"), Comment("Minimum radial distance between calo cluster and helix."), 0.};
      fhicl::Atom<double> MinDeltaT {
        Name("MinDeltaT"), Comment("Minimum time between calo cluster and helix."), 0.};
      fhicl::Atom<double> MinCR {
        Name("MinCR"), Comment("Minimum radial coordinate of the cluster."), 0.};
      fhicl::Atom<double> MaxCT {
        Name("MaxCT"), Comment("Maximum time of the cluster."), 2000.};
      fhicl::Atom<double> MinCT {
        Name("MinCT"), Comment("Minimum time of the cluster."), 0.};
      fhicl::Atom<int> MaxNCrys {
        Name("MaxNCrys"), Comment("Maximum number of crystals in a cluster."), 100};
      fhicl::Atom<double> MinCalE {
        Name("MinCalE"), Comment("Minimum Calo Cluster Energy."), 0};
      fhicl::Atom<double> MinCalERatio {
        Name("MinCalERatio"), Comment("Minimum Primary Hit to Total Cluster Energy Ratio."), 0.0};
      fhicl::Atom<double> MinCalE2Ratio {
        Name("MinCalE2Ratio"), Comment("Minimum 2 Primary Hits to Total Cluster Energy Ratio."), 0.0};
      fhicl::Atom<double> MinCalE9Ratio {
        Name("MinCalE9Ratio"), Comment("Minimum Ratio of energy of 9 crystals centered on primary to  Total Cluster Energy."), 0.0};
      fhicl::Atom<double> MinCalE25Ratio {
        Name("MinCalE25Ratio"), Comment("Minimum Ratio of energy of 25 crystals centered on primary to  Total Cluster Energy."), 0.0};
      fhicl::Atom<double> MinCHERatio {
        Name("MinCHERatio"), Comment("Minimum Calo Helix Energy Ratio."), 0.0};
      fhicl::Atom<double> MinHelP {
        Name("MinHelP"), Comment("Minimum Helix Momentum"), 0.0};
      fhicl::Atom<double> MaxHelChi2dXY {
        Name("MaxHelChi2dXY"), Comment("Maximum helix Chi2 in dXY"), 100.0};
      fhicl::Atom<double> MaxHelChi2dZPhi {
        Name("MaxHelChi2dZPhi"), Comment("Maximum helix Chi2 in dZPhi"), 100.0};
      fhicl::Atom<bool> ApplyCuts {
        Name("ApplyCuts"), Comment("Boolean Variable for whether to apply cuts."), false};


    };

    using Parameters = art::EDFilter::Table<Config>;
    explicit PhotonFilter(const Parameters& conf);

  private:

    // Histogram struct
    struct hist_c {
      TH1F* _hCalE;
      TH1F* _hSimCalE;
      TH1F* _hSimE;
      TH1F* _hCalERatio;

      TH1F* _hCalECrys;
      TH1F* _hZRes_CloseHelix;
      TH1F* _hCRes_CloseHelix;
      TH1F* _hTRes;
      TH1F* _hNCrystals;

      TH1F* _hZcoord_calo;
      TH1F* _hHelixP;
      TH1F* _hHelixChi2_dXY;
      TH1F* _hHelixChi2_dZPhi;
      TH1F* _hHelixRadius;
      TH1F* _hClusterR;
      TH1F* _hClusterT;
      TH1F* _hSimPdg;
      TH1F* _hClusterHelixEnergyRatio;
      TH1F* _hCalE2Ratio;
      TH1F* _hCalE9Ratio;
      TH1F* _hCalE25Ratio;
      TH1F* _hClusterSecondMoment;
    };
    struct hist_e {
      TH1F* _hNClusters;
    };
    enum {kNhists = 100};
    hist_c cluster_hists [kNhists];
    hist_e event_hists [kNhists];

    bool filter(art::Event& event) override;
    bool reportResiduals(const CaloCluster& cluster, const HelixSeedCollection& helices, const CaloClusterMC& caloClusterTruth);
    bool beginRun(art::Run& run) override;
    void bookClusterHistograms(hist_c& cluster_hists, art::TFileDirectory& Dir);
    void bookEventHistograms(hist_e& event_hists, art::TFileDirectory& Dir);



    art::InputTag _ccTag;
    art::InputTag _hsTag;
    art::InputTag _ccMCTag;
    art::InputTag _ppTag;
    int _debugLevel;
    double _mindR;
    double _mindT;
    double _minCR;
    double _maxCT;
    double _minCT;
    int _maxNCrys;
    double _minCalE;
    double _minCalERatio;
    double _minCalE2Ratio;
    double _minCalE9Ratio;
    double _minCalE25Ratio;

    double _minCHERatio;
    double _minHelP;
    double _maxHelChi2dXY;
    double _maxHelChi2dZPhi;
    bool _applyCuts;

    const mu2e::Calorimeter*             _calorimeter;
    const mu2e::SimParticle*             _primaryParticle;


  };

  PhotonFilter::PhotonFilter(const Parameters& conf)
    : art::EDFilter{conf}
    , _ccTag(conf().CaloClusterCollection())
    , _hsTag(conf().HelixSeedCollection())
    , _ccMCTag(conf().CaloClusterMC())
    , _ppTag(conf().PrimaryParticle())
    , _debugLevel(conf().DebugLevel())
    , _mindR(conf().MinDeltaR())
    , _mindT(conf().MinDeltaT())
    , _minCR(conf().MinCR())
    , _maxCT(conf().MaxCT())
    , _minCT(conf().MinCT())
    , _maxNCrys(conf().MaxNCrys())
    , _minCalE(conf().MinCalE())
    , _minCalERatio(conf().MinCalERatio())
    , _minCalE2Ratio(conf().MinCalE2Ratio())
    , _minCalE9Ratio(conf().MinCalE9Ratio())
    , _minCalE25Ratio(conf().MinCalE25Ratio())
    , _minCHERatio(conf().MinCHERatio())
    , _minHelP(conf().MinHelP())
    , _maxHelChi2dXY(conf().MaxHelChi2dXY())
    , _maxHelChi2dZPhi(conf().MaxHelChi2dZPhi())
    , _applyCuts(conf().ApplyCuts())



  {
    produces<TriggerInfo>();
    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory DirE0 =  tfs->mkdir("all_events");
    art::TFileDirectory DirE1 =  tfs->mkdir("passed_events");

    bookEventHistograms(event_hists[0], DirE0);
    bookEventHistograms(event_hists[1], DirE1);

    art::TFileDirectory Dir0 =  tfs->mkdir("all_clusters");
    bookClusterHistograms(cluster_hists[0], Dir0);
    art::TFileDirectory Dir1 =  tfs->mkdir("all_photon_clusters");
    bookClusterHistograms(cluster_hists[1], Dir1);
    art::TFileDirectory Dir2 =  tfs->mkdir("all_background_clusters");
    bookClusterHistograms(cluster_hists[2], Dir2);
    art::TFileDirectory Dir10 =  tfs->mkdir("C_passed_clusters");
    bookClusterHistograms(cluster_hists[10], Dir10);
    art::TFileDirectory Dir20 =  tfs->mkdir("C_rejected_clusters");
    bookClusterHistograms(cluster_hists[20], Dir20);
    art::TFileDirectory Dir11 =  tfs->mkdir("C_passed_photon_clusters");
    bookClusterHistograms(cluster_hists[11], Dir11);
    art::TFileDirectory Dir21 =  tfs->mkdir("C_rejected_photon_clusters");
    bookClusterHistograms(cluster_hists[21], Dir21);
    art::TFileDirectory Dir12 =  tfs->mkdir("C_passed_background_clusters");
    bookClusterHistograms(cluster_hists[12], Dir12);
    art::TFileDirectory Dir22 =  tfs->mkdir("C_rejected_background_clusters");
    bookClusterHistograms(cluster_hists[22], Dir22);

  }

  void PhotonFilter::bookEventHistograms(hist_e& hists, art::TFileDirectory& Dir){
    TH1::AddDirectory(0);
    hists._hNClusters = Dir.make<TH1F>("hNClusters",   "Number of Clusters", 100, 0, 100);
  }
  void PhotonFilter::bookClusterHistograms(hist_c& hists, art::TFileDirectory& Dir){
    TH1::AddDirectory(0);
    hists._hCalE   = Dir.make<TH1F>("hCalE",   "Cluster energy, MeV/c", 250, 0, 250);
    hists._hSimCalE   = Dir.make<TH1F>("hSimCalE",   "Simulated Cluster energy, MeV/c", 125, 0, 250);
    hists._hSimE   = Dir.make<TH1F>("hSimE",   "Simulated particle energy, MeV/c", 125, 0, 250);

    hists._hCalERatio   = Dir.make<TH1F>("hCalERatio",   "Ratio of energy in the primary crystal vs total energy, MeV/c", 100, 0, 1);
    hists._hCalE2Ratio   = Dir.make<TH1F>("hCalE2Ratio",   "Ratio of energy in the two highest energy crystals vs total energy, MeV/c", 100, 0, 1);
    hists._hCalE9Ratio   = Dir.make<TH1F>("hCalE9Ratio",   "Ratio of energy of 9 crystals centered on  highest energy crystal vs total energy, MeV/c", 100, 0, 1);
    hists._hCalE25Ratio   = Dir.make<TH1F>("hCalE25Ratio",   "Ratio of energy of 25 crystals centered on  highest energy crystal vs total energy, MeV/c", 100, 0, 1);

    hists._hCalECrys   = Dir.make<TH1F>("hCalECrys",   "Cluster energy per Crystal, MeV/c", 125, 0, 250);

    hists._hCRes_CloseHelix   = Dir.make<TH1F>("hCResidual_CloesestHelix",   "Helix of lowest Residual, extrapolated from Center Coordinate of Cluster (mm)", 500, 0, 1000);
    hists._hTRes   = Dir.make<TH1F>("hTRes",   "Difference in time between helix and cluster (ns)", 1000, 0, 2000);

    hists._hNCrystals   = Dir.make<TH1F>("hNCrystals",   "Number of crystals in a cluster", 100, 0, 10);

    hists._hHelixP = Dir.make<TH1F>("hHelixP",  "Helix Momentum (MeV)",     250, 0, 500);
    hists._hHelixChi2_dXY = Dir.make<TH1F>("hHelixChi2_dXY",  "Helix Chi2 in dXY",     250, 0, 50);
    hists._hHelixChi2_dZPhi = Dir.make<TH1F>("hHelixChi2_dZPhi",  "Helix Chi2 in dZPhi",     100, 0, 50);
    hists._hHelixRadius = Dir.make<TH1F>("hHelixRadius",  "Helix Radius",     250, 0, 500);

    hists._hClusterR = Dir.make<TH1F>("clusterR",  "Cluster Radial Coordinate (mm)",     500, 0, 1000);
    hists._hClusterT = Dir.make<TH1F>("clusterT",  "Cluster Time (ns)",     1000, 0, 2000);

    hists._hSimPdg = Dir.make<TH1F>("SimPdgID",  "Sim PDG ID",     1250, 0, 2500);
    hists._hClusterHelixEnergyRatio = Dir.make<TH1F>("ClusterHelixEnergyRatio",  "Ratio of Helix to Cluster Energy",     100, 0, 1);
    hists._hClusterSecondMoment = Dir.make<TH1F>("ClusterSecondMoment",  "Cluster Second Moment",     200, 0, 20000);





  }

  bool PhotonFilter::beginRun(art::Run& run) {
    GeomHandle<mu2e::Calorimeter> ch;
    _calorimeter = ch.get();
    if (_debugLevel >= 1) {
      std::cout << "[PhotonFilter] beginRun(): Calorimeter geometry loaded." << std::endl;
    }
    return true;
  }


  bool PhotonFilter::filter(art::Event& evt) {

    auto hsH = evt.getValidHandle<HelixSeedCollection>(_hsTag);
    auto ccH = evt.getValidHandle<CaloClusterCollection>(_ccTag);
    auto ccMC = evt.getValidHandle<CaloClusterMCCollection>(_ccMCTag);
    auto pp = evt.getValidHandle<PrimaryParticle>(_ppTag);

    const auto& helices = *hsH;
    const auto& clusters = *ccH;
    const auto& caloClusterTruth = *ccMC;
    const auto& primaryParticles = *pp;
    const auto& simps = primaryParticles.primarySimParticles();

    _primaryParticle = nullptr;
    if (!simps.empty() && simps.front().isNonnull()) {
      _primaryParticle = simps.front().get();   // safe: returns nullptr if somehow null, but isNonnull already checked
    }
    auto triginfo = std::make_unique<TriggerInfo>();


    event_hists[0]._hNClusters->Fill(clusters.size());


    // Add all clusters as art::Ptr to TriggerInfo
    for (size_t i = 0; i < clusters.size(); ++i) {
      if (reportResiduals(clusters.at(i), helices, caloClusterTruth[i])){
        art::Ptr<CaloCluster> clusterPtr(ccH, i);
        triginfo->_caloClusters.push_back(clusterPtr);
      }
    }
    const bool passed = !triginfo->_caloClusters.empty();
    if (passed){
      event_hists[1]._hNClusters->Fill(triginfo->_caloClusters.size());}

    evt.put(std::move(triginfo));
    return passed ;
  }




  bool PhotonFilter::reportResiduals(const CaloCluster& cluster, const HelixSeedCollection& helices, const CaloClusterMC& caloClusterTruth) {


    CLHEP::Hep3Vector gpos = _calorimeter->geomUtil().diskToMu2e(cluster.diskID(), cluster.cog3Vector());

    CLHEP::Hep3Vector tpos = _calorimeter->geomUtil().mu2eToTracker(gpos);

    double offset = _calorimeter->caloInfo().getDouble("diskCaseZLength");
    offset += _calorimeter->caloInfo().getDouble("BPPipeZOffset");
    offset += _calorimeter->caloInfo().getDouble("BPHoleZLength");
    offset += _calorimeter->caloInfo().getDouble("FEEZLength");
    offset /= 2.0;

    XYZVectorF clusterPos(tpos.x(), tpos.y(), tpos.z() - offset);

    const auto clusterEDep = cluster.energyDep();

    ClusterUtils clu(*_calorimeter, cluster);

    const double e2  = clu.e2()/clusterEDep;
    const double e9  = clu.e9()/clusterEDep;
    const double e25 = clu.e25()/clusterEDep;
    const double secMoment = clu.secondMoment();

    const int  clusterHits          = cluster.caloHitsPtrVector().size();
    bool bestHelixExists = false;


    float cR = sqrt(clusterPos.x()*clusterPos.x()+clusterPos.y()*clusterPos.y());
    float cT = cluster.time();

    const auto eDepMCs =  caloClusterTruth.energyDeposits();

    std::unordered_map<int,double> eDepSimMap;
    // find PDG with largest total EDep
    double maxEDep = -1.0;
    art::Ptr<SimParticle> mainSim;


    // accumulate energy by sim ID and PDG
    for (size_t j = 0; j < eDepMCs.size(); ++j) {
      const auto& eDepMC = eDepMCs[j];
      art::Ptr<SimParticle> sim = eDepMC.sim();
      const int simID = sim->id().asInt();
      eDepSimMap[simID] += eDepMC.energyDep();
      if (eDepSimMap[simID] > maxEDep) {
        maxEDep =  eDepSimMap[simID];
        mainSim = sim;

      }
    }

    int mainPdg = 0;
    bool isPrimaryRel = false;

    if(mainSim.isNull()) {
      std::cout << "[PhotonFilter::" << __func__ << "] No main sim found! N(edeps) = " << eDepMCs.size() << std::endl;
    } else {

      mainPdg = mainSim->pdgId();

      if(_primaryParticle != nullptr) {
        art::Ptr<SimParticle> sim_parent = mainSim;
        while(sim_parent.isNonnull()) {
          if(&(*_primaryParticle) == &(*sim_parent)) {
            isPrimaryRel = true;
            break;
          }
          sim_parent = sim_parent.get()->parent();
        }
      }
    }

    double clusterSimE = maxEDep;
    double genE0 = 0;
    if(_primaryParticle != nullptr) {genE0 =  _primaryParticle->startMomentum().e();}

    float minCResidual = std::numeric_limits<float>::max();
    double bestHelixMom = 0;
    double bestHelixChi2dXY = 0;
    double bestHelixChi2dZPhi = 0;
    double bestHelixRadius = 0;
    double bestHelixT = -1;


    for (size_t i = 0; i < helices.size(); ++i) {
      const auto& helix = helices[i];
      const auto& hval = helix.helix();

      if (hval.momentum()*0.3 < 90) continue;

      // Extract helix center and radius
      float radius = hval.radius();
      float centerX = hval.centerx();
      float centerY = hval.centery();

      float deltaCR = std::fabs(std::hypot(clusterPos.x() - centerX,
                                           clusterPos.y() - centerY) - radius);
      if (deltaCR < minCResidual) {
        minCResidual = deltaCR;
        bestHelixExists = true;
        bestHelixT = helix.t0()._t0;
        bestHelixMom = hval.momentum()*0.3;
        bestHelixRadius = hval.radius();
        bestHelixChi2dXY = hval.chi2dXY();
        bestHelixChi2dZPhi = hval.chi2dZPhi();
      }
    }

    double bestdT = std::abs(bestHelixT - cT);
    int index = 0;

    cluster_hists[index]._hCalE->Fill(clusterEDep);
    cluster_hists[index]._hSimCalE->Fill(clusterSimE);
    cluster_hists[index]._hSimE->Fill(genE0);
    cluster_hists[index]._hNCrystals->Fill(clusterHits);
    cluster_hists[index]._hCalECrys->Fill(clusterEDep/clusterHits);
    cluster_hists[index]._hCalERatio->Fill(cluster.caloHitsPtrVector()[0]->energyDep()/clusterEDep);
    cluster_hists[index]._hCalE2Ratio->Fill(e2);
    cluster_hists[index]._hCalE9Ratio->Fill(e9);
    cluster_hists[index]._hCalE25Ratio->Fill(e25);
    cluster_hists[index]._hClusterSecondMoment->Fill(secMoment);



    cluster_hists[index]._hClusterR->Fill(cR);
    cluster_hists[index]._hClusterT->Fill(cT);
    cluster_hists[index]._hSimPdg->Fill(mainPdg);


    if (bestHelixExists) {
      cluster_hists[index]._hCRes_CloseHelix->Fill(minCResidual);
      cluster_hists[index]._hTRes->Fill(bestdT);
      cluster_hists[index]._hHelixP->Fill(bestHelixMom);
      cluster_hists[index]._hHelixRadius->Fill(bestHelixRadius);
      cluster_hists[index]._hHelixChi2_dXY->Fill(bestHelixChi2dXY);
      cluster_hists[index]._hHelixChi2_dZPhi->Fill(bestHelixChi2dZPhi);
      cluster_hists[index]._hClusterHelixEnergyRatio->Fill(clusterEDep/bestHelixMom);

    }

    if (isPrimaryRel) {
      index = 1;
    }
    else {index = 2;}

    cluster_hists[index]._hCalE->Fill(clusterEDep);
    cluster_hists[index]._hSimE->Fill(genE0);
    cluster_hists[index]._hSimCalE->Fill(clusterSimE);
    cluster_hists[index]._hNCrystals->Fill(clusterHits);
    cluster_hists[index]._hCalECrys->Fill(clusterEDep / clusterHits);
    cluster_hists[index]._hCalERatio->Fill(cluster.caloHitsPtrVector()[0]->energyDep() / clusterEDep);
    cluster_hists[index]._hCalE2Ratio->Fill(e2);
    cluster_hists[index]._hCalE9Ratio->Fill(e9);
    cluster_hists[index]._hCalE25Ratio->Fill(e25);
    cluster_hists[index]._hClusterSecondMoment->Fill(secMoment);


    cluster_hists[index]._hClusterR->Fill(cR);
    cluster_hists[index]._hClusterT->Fill(cT);

    cluster_hists[index]._hSimPdg->Fill(mainPdg);

    if (bestHelixExists) {
      cluster_hists[index]._hCRes_CloseHelix->Fill(minCResidual);
      cluster_hists[index]._hTRes->Fill(bestdT);
      cluster_hists[index]._hHelixP->Fill(bestHelixMom);
      cluster_hists[index]._hHelixRadius->Fill(bestHelixRadius);
      cluster_hists[index]._hHelixChi2_dXY->Fill(bestHelixChi2dXY);
      cluster_hists[index]._hHelixChi2_dZPhi->Fill(bestHelixChi2dZPhi);
      cluster_hists[index]._hClusterHelixEnergyRatio->Fill(clusterEDep/bestHelixMom);

    }


    bool accepted = true;
    if (_applyCuts){
      accepted = (
                  clusterHits   < _maxNCrys &&
                  clusterEDep   > _minCalE &&
                  cR > _minCR && cT > _minCT && cT < _maxCT &&
                  cluster.caloHitsPtrVector()[0]->energyDep() / clusterEDep > _minCalERatio && e2 > _minCalE2Ratio && e9 > _minCalE9Ratio && e25 > _minCalE25Ratio                    );
      if (bestHelixExists) {accepted = (accepted && minCResidual > _mindR && bestdT > _mindT && bestHelixMom > _minHelP && bestHelixChi2dXY < _maxHelChi2dXY && bestHelixChi2dZPhi < _maxHelChi2dZPhi && clusterEDep/bestHelixMom > _minCHERatio);}
    }
    int indexAccept = 0;
    int indexPdg = 0;
    if (accepted) {
      indexAccept = 10;
      if (isPrimaryRel) {indexPdg = 11;}
      else {indexPdg = 12;}
    }
    else{indexAccept = 20;
      if (isPrimaryRel) {indexPdg = 21;}
      else {indexPdg = 22;}
    }
    if (indexAccept > 0) {
      cluster_hists[indexAccept]._hCalE->Fill(clusterEDep);
      cluster_hists[indexAccept]._hSimCalE->Fill(clusterSimE);
      cluster_hists[indexAccept]._hSimE->Fill(genE0);
      cluster_hists[indexAccept]._hNCrystals->Fill(clusterHits);
      cluster_hists[indexAccept]._hCalECrys->Fill(clusterEDep/clusterHits);
      cluster_hists[indexAccept]._hCalERatio->Fill(cluster.caloHitsPtrVector()[0]->energyDep()/clusterEDep);
      cluster_hists[indexAccept]._hCalE2Ratio->Fill(e2);
      cluster_hists[indexAccept]._hCalE9Ratio->Fill(e9);
      cluster_hists[indexAccept]._hCalE25Ratio->Fill(e25);
      cluster_hists[indexAccept]._hClusterSecondMoment->Fill(secMoment);


      cluster_hists[indexAccept]._hClusterR->Fill(cR);
      cluster_hists[indexAccept]._hClusterT->Fill(cT);
      cluster_hists[indexAccept]._hSimPdg->Fill(mainPdg);

      if (bestHelixExists) {
        cluster_hists[indexAccept]._hCRes_CloseHelix->Fill(minCResidual);
        cluster_hists[indexAccept]._hTRes->Fill(bestdT);
        cluster_hists[indexAccept]._hHelixP->Fill(bestHelixMom);
        cluster_hists[indexAccept]._hHelixRadius->Fill(bestHelixRadius);

        cluster_hists[indexAccept]._hHelixChi2_dXY->Fill(bestHelixChi2dXY);
        cluster_hists[indexAccept]._hHelixChi2_dZPhi->Fill(bestHelixChi2dZPhi);
        cluster_hists[indexAccept]._hClusterHelixEnergyRatio->Fill(clusterEDep/bestHelixMom);

      }
    }

    if (indexPdg > 0) {
      cluster_hists[indexPdg]._hCalE->Fill(clusterEDep);
      cluster_hists[indexPdg]._hSimCalE->Fill(clusterSimE);
      cluster_hists[indexPdg]._hSimE->Fill(genE0);
      cluster_hists[indexPdg]._hNCrystals->Fill(clusterHits);
      cluster_hists[indexPdg]._hCalECrys->Fill(clusterEDep / clusterHits);
      cluster_hists[indexPdg]._hCalERatio->Fill(cluster.caloHitsPtrVector()[0]->energyDep() / clusterEDep);
      cluster_hists[indexPdg]._hCalE2Ratio->Fill(e2);
      cluster_hists[indexPdg]._hCalE9Ratio->Fill(e9);
      cluster_hists[indexPdg]._hCalE25Ratio->Fill(e25);
      cluster_hists[indexPdg]._hClusterSecondMoment->Fill(secMoment);


      cluster_hists[indexPdg]._hClusterR->Fill(cR);
      cluster_hists[indexPdg]._hClusterT->Fill(cT);
      cluster_hists[indexPdg]._hSimPdg->Fill(mainPdg);

      if (bestHelixExists) {
        cluster_hists[indexPdg]._hCRes_CloseHelix->Fill(minCResidual);
        cluster_hists[indexPdg]._hTRes->Fill(bestdT);
        cluster_hists[indexPdg]._hHelixP->Fill(bestHelixMom);
        cluster_hists[indexPdg]._hHelixRadius->Fill(bestHelixRadius);
        cluster_hists[indexPdg]._hHelixChi2_dXY->Fill(bestHelixChi2dXY);
        cluster_hists[indexPdg]._hHelixChi2_dZPhi->Fill(bestHelixChi2dZPhi);
        cluster_hists[indexPdg]._hClusterHelixEnergyRatio->Fill(clusterEDep/bestHelixMom);

      }
    }

    if (_debugLevel >= 1) {
      std::cout << " Photon Filter: best match c residual = " << minCResidual << std::endl;

    }
    if (_applyCuts){
      return accepted;}
    else {return true;}

  }
}
// namespace mu2e

using mu2e::PhotonFilter;
DEFINE_ART_MODULE(PhotonFilter)
