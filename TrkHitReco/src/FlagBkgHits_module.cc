#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

#include "Offline/ConditionsService/inc/ConditionsHandle.hh"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"

#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/BkgCluster.hh"
#include "Offline/RecoDataProducts/inc/BkgClusterHit.hh"

#include "Offline/TrkHitReco/inc/TNTClusterer.hh"
#include "Offline/TrkHitReco/inc/ScanClusterer.hh"
#include "Offline/TrkHitReco/inc/TrainBkgDiag.hxx"

#include <string>
#include <vector>


//NEW CLASS
namespace TMVA_SOFIE_TrainBkgDiag {
  class Session;
}

namespace mu2e
{

  class FlagBkgHits : public art::EDProducer
  {
    public:

      struct Config
      {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;
        fhicl::Atom<art::InputTag>            comboHitCollection{   Name("ComboHitCollection"),   Comment("ComboHit collection name") };
        fhicl::Atom<art::InputTag>            strawHitCollection{   Name("StrawHitCollection"),   Comment("StrawHit collection name") };
        fhicl::Atom<unsigned>                 minActiveHits{        Name("MinActiveHits"),        Comment("Minumim number of active hits in a cluster") };
        fhicl::Atom<unsigned>                 minNPlanes{           Name("MinNPlanes"),           Comment("Minumim number of planes in a cluster") };
        fhicl::Atom<float>                    clusterPositionError{ Name("ClusterPositionError"), Comment("Cluster poisiton error") };
        fhicl::Atom<int>                      clusterAlgorithm{     Name("ClusterAlgorithm"),     Comment("Clusterer algorithm") };
        fhicl::Atom<bool>                     filterHits{           Name("FilterHits"),           Comment("Produce filtered ComboHit collection")  };
        fhicl::Atom<bool>                     flagComboHits {       Name("FlagComboHits"),        Comment("Produce ComboHit-levelflag collection") };
        fhicl::Atom<bool>                     flagStrawHits {       Name("FlagStrawHits"),        Comment("Produce StrawHit-level flag collection") };
        fhicl::Sequence<std::string>          backgroundMask{       Name("BackgroundMask"),       Comment("Bkg hit selection mask") };
        fhicl::Sequence<std::string>          stereoSelection{      Name("StereoSelection"),      Comment("Stereo hit selection mask") };
        fhicl::Atom<bool>                     saveBkgClusters{      Name("SaveBkgClusters"),      Comment("Save bkg clusters") };
        fhicl::Atom<int>                      debugLevel{           Name("DebugLevel"),           Comment("Debug"),0 };
        fhicl::Table<TNTClusterer::Config>    TNTClustering{        Name("TNTClustering"),        Comment("TNT Clusterer config") };
        fhicl::Table<ScanClusterer::Config>   ScanClustering{       Name("ScanClustering"),       Comment("Scan Clusterer config") };
        fhicl::Atom<float>                    kerasQuality{         Name("KerasQuality"),         Comment("Keras quality cut") };
      };

      enum clusterer {TwoNiveauThreshold=1, ComptonKiller=2};
      explicit FlagBkgHits(const art::EDProducer::Table<Config>& config);
      void beginJob() override;
      void produce(art::Event& event) override;

    private:
      const art::ProductToken<ComboHitCollection> chtoken_;
      const art::ProductToken<StrawHitCollection> shtoken_;
      unsigned                                    minnhits_;
      unsigned                                    minnp_;
      bool                                        filter_, flagsh_, flagch_;
      bool                                        savebkg_;
      StrawHitFlag                                bkgmsk_, stereo_;
      BkgClusterer*                               clusterer_;
      float                                       cperr2_;
      int const                                   debug_;
      float                                       kerasQ_;
      int                                         iev_;
      std::shared_ptr<TMVA_SOFIE_TrainBkgDiag::Session> sofiePtr;

      void classifyCluster(BkgClusterCollection& bkgccolFast, BkgClusterCollection& bkgccol,
          StrawHitFlagCollection& chfcol, const ComboHitCollection& chcol) const;
      void countHits(      const BkgCluster& cluster, unsigned& nactive, unsigned& nstereo, const ComboHitCollection& chcol) const;
      void countPlanes(    const BkgCluster& cluster, std::vector<float>& kerasvars, const ComboHitCollection& chcol) const;
      int  findClusterIdx( BkgClusterCollection& bkgccol, unsigned ich) const;
  };


  FlagBkgHits::FlagBkgHits(const art::EDProducer::Table<Config>& config) :
    art::EDProducer{config},
    chtoken_{     consumes<ComboHitCollection>(config().comboHitCollection()) },
    shtoken_{     consumes<StrawHitCollection>(config().strawHitCollection()) },
    minnhits_(    config().minActiveHits() ),
    minnp_(       config().minNPlanes()),
    filter_(      config().filterHits()),
    flagsh_(      config().flagStrawHits()),
    flagch_(      config().flagComboHits()),
    savebkg_(     config().saveBkgClusters()),
    bkgmsk_(      config().backgroundMask()),
    stereo_(      config().stereoSelection()),
    debug_(       config().debugLevel()),
    kerasQ_(      config().kerasQuality()),
    iev_(0)
    {
      // Must call consumesMany because fillStrawHitIndices calls getManyByType.
      consumesMany<ComboHitCollection>();
      sofiePtr = std::make_shared<TMVA_SOFIE_TrainBkgDiag::Session>("Offline/TrkHitReco/data/TrainBkgDiag.dat");

      if (flagsh_) produces<StrawHitFlagCollection>("StrawHits");
      if (flagch_) produces<StrawHitFlagCollection>("ComboHits");
      if (filter_) produces<ComboHitCollection>();

      if (savebkg_)
      {
        produces<BkgClusterHitCollection>();
        produces<BkgClusterCollection>();
      }
      float cperr = config().clusterPositionError();
      cperr2_ = cperr*cperr;

      clusterer ctype = static_cast<clusterer>(config().clusterAlgorithm());
      switch ( ctype )
      {
        case TwoNiveauThreshold:
          clusterer_ = new TNTClusterer(config().TNTClustering());
          break;
        case ComptonKiller:
          clusterer_ = new ScanClusterer(config().ScanClustering());
          break;
        default:
          throw cet::exception("RECO")<< "Unknown clusterer" << ctype << std::endl;
      }
    }


  void FlagBkgHits::beginJob()
  {
    clusterer_->init();
  }


  //------------------------------------------------------------------------------------------
  void FlagBkgHits::produce(art::Event& event )
  {

    auto chH = event.getValidHandle(chtoken_);
    const ComboHitCollection& chcol = *chH.product();
    unsigned nch = chcol.size();


    // the primary output is either a deep copy of selected inputs or a flag collection on those
    // intermediate results: keep these on the heap unless requested for diagnostics later
    BkgClusterCollection bkgccol,bkgccolFast;
    BkgClusterHitCollection bkghitcol;
    bkgccol.reserve(nch/2);
    if (savebkg_) bkghitcol.reserve(nch);


    // find clusters, sort is needed for recovery algorithm. bkgccolFast has hits that are autmoatically marked as bkg.
    clusterer_->findClusters(bkgccolFast,bkgccol,chcol, iev_);
    sort(bkgccol.begin(),bkgccol.end(),[](const BkgCluster& c1,const BkgCluster& c2) {return c1.time() < c2.time();});


    // classify clusters
    StrawHitFlagCollection chfcol(nch);
    //classifyCluster(bkgccolFast, bkgccol, bkgqcol, chfcol,chcol);
    classifyCluster(bkgccolFast, bkgccol, chfcol, chcol);

    //produce BkgClusterHit info collection
    if (savebkg_)
    {
      for (size_t ich=0;ich < chcol.size(); ++ich)
      {
        const ComboHit& ch = chcol[ich];
        int icl = findClusterIdx(bkgccol,ich);
        if (icl > -1) bkghitcol.emplace_back(BkgClusterHit(clusterer_->distance(bkgccol[icl],ch),ch.flag()));
        else          bkghitcol.emplace_back(BkgClusterHit(999.0,ch.flag()));
      }
    }

    //produce filtered ComboHit collection
    if (filter_)
    {
      auto chcolFilter = std::make_unique<ComboHitCollection>();
      chcolFilter->reserve(nch);
      // same parent as the original collection
      chcolFilter->setParent(chcol.parent());
      for(size_t ich=0;ich < nch; ++ich)
      {
        StrawHitFlag const& flag = chfcol[ich];
        if (!flag.hasAnyProperty(bkgmsk_))
        {
          chcolFilter->push_back(chcol[ich]);
          chcolFilter->back()._flag.merge(flag);
        }
      }
      event.put(std::move(chcolFilter));
    }

    //produce StrawHit flags
    // Note: fillStrawHitIndices is quite slow and should be improved
    if (flagsh_)
    {
      auto shH = event.getValidHandle(shtoken_);
      const StrawHitCollection* shcol  = shH.product();

      unsigned nsh = shcol->size();
      auto shfcol  = std::make_unique<StrawHitFlagCollection>(nsh);
      std::vector<std::vector<StrawHitIndex> > shids;
      chcol.fillStrawHitIndices(event,shids);
      for (size_t ich = 0;ich < nch;++ich)
      {
        StrawHitFlag flag = chfcol[ich];
        flag.merge(chcol[ich].flag());
        for(auto ish : shids[ich]) (*shfcol)[ish] = flag;
      }

      event.put(std::move(shfcol),"StrawHits");
    }

    //produce ComboHit flags
    if (flagch_)
    {
      for(size_t ich=0;ich < nch; ++ich) chfcol[ich].merge(chcol[ich].flag());
      event.put(std::make_unique<StrawHitFlagCollection>(std::move(chfcol)),"ComboHits");
    }

    //produce background collection
    if (savebkg_)
    {
      event.put(std::make_unique<BkgClusterHitCollection>(bkghitcol));
      event.put(std::make_unique<BkgClusterCollection>(bkgccol));
    }

    ++iev_;
    return;
  }


  //------------------------------------------------------------------------------------------
  void FlagBkgHits::classifyCluster(BkgClusterCollection& bkgccolFast, BkgClusterCollection& bkgccol,
      StrawHitFlagCollection& chfcol, const ComboHitCollection& chcol) const
  {

    for (auto& cluster : bkgccolFast)
    {
      StrawHitFlag flag(StrawHitFlag::bkgclust);
      flag.merge(StrawHitFlag::bkg);
      for (const auto& chit : cluster.hits()) chfcol[chit] = flag;
    }

    for (auto& cluster : bkgccol)
    {

      unsigned nactive, nstereo;
      countHits(cluster, nactive, nstereo,chcol);
      unsigned nhits = nstereo > 0 ? nstereo : nactive;
      if (nhits < minnhits_) return;
      std::vector<float> kerasvars(11,0.0);
      kerasvars.at(0) = sqrtf(cluster.pos().perp2());
      kerasvars.at(6) = nhits;
      countPlanes(cluster,kerasvars,chcol);

      if (kerasvars.at(4) >= minnp_)
      {
        std::vector<float> hz;
        for (const auto& chit : cluster.hits()) hz.push_back(chcol[chit].pos().z());

        // find the min, max and largest gap from the sorted Z positions
        std::sort(hz.begin(),hz.end());
        float zgap = 0.0;
        for (unsigned iz=1;iz<hz.size();++iz) zgap=std::max(zgap,hz[iz]-hz[iz-1]);
        kerasvars.at(1) = hz.front(); // Z min
        kerasvars.at(2) = hz.back(); // Z max
        kerasvars.at(3) = zgap; // max Z gap
      }
      else
      {
        kerasvars.at(1) = -1.0;
        kerasvars.at(2)  = -1.0;
        kerasvars.at(3)  = -1.0;
      }

      double sumEdep(0.);
      double sqrSumDeltaTime(0.);
      double sqrSumDeltaX(0.);
      double sqrSumDeltaY(0.);
      for (const auto& chit : cluster.hits())
      {
        sumEdep +=  chcol[chit].energyDep()/chcol[chit].nStrawHits();
        sqrSumDeltaX += std::pow(chcol[chit].pos().x() - cluster.pos().x(),2);
        sqrSumDeltaY += std::pow(chcol[chit].pos().y() - cluster.pos().y(),2);
        sqrSumDeltaTime += std::pow(chcol[chit].time() - cluster.time(),2);
      }

      kerasvars.at(7) = sumEdep/nhits;
      kerasvars.at(8) = std::sqrt(sqrSumDeltaX/nhits);
      kerasvars.at(9) = std::sqrt(sqrSumDeltaY/nhits);
      kerasvars.at(10) = std::sqrt(sqrSumDeltaTime/nhits);

      auto kerasout = sofiePtr->infer(kerasvars.data());

      StrawHitFlag flag(StrawHitFlag::bkgclust);
      if (kerasout[0] > kerasQ_)
      {
        flag.merge(StrawHitFlag::bkg);
        if (savebkg_)
        {
          cluster._flag.merge(BkgClusterFlag::bkg);
          if(nstereo > 0) cluster._flag.merge(BkgClusterFlag::stereo);
        }
      }

      for (const auto& chit : cluster.hits()) chfcol[chit] = flag;
    }
  }

  //-------------------------------------------------------------------------------------------------------------
  void FlagBkgHits::countPlanes(const BkgCluster& cluster, std::vector<float> & kerasvars, const ComboHitCollection& chcol) const
  {
    std::array<int,StrawId::_nplanes> hitplanes{0};
    for (const auto& chit : cluster.hits())
    {
      const ComboHit& ch = chcol[chit];
      hitplanes[ch.strawId().plane()] += ch.nStrawHits();
    }

    unsigned ipmin(0),ipmax(StrawId::_nplanes-1);
    while (hitplanes[ipmin]==0) ++ipmin;
    while (hitplanes[ipmax]==0) --ipmax;

    unsigned npexp(0),np(0),nphits(0);
    for(unsigned ip = ipmin; ip <= ipmax; ++ip)
    {
      npexp++; // should use TTracker to see if plane is physically present FIXME!
      if (hitplanes[ip]> 0)++np;
      nphits += hitplanes[ip];
    }

    kerasvars.at(4) = np;// # of planes
    kerasvars.at(5) = static_cast<float>(np)/static_cast<float>(npexp);// fraction of planes

  }


  //----------------------------------------------------------------------------------------------------------------------------------
  void FlagBkgHits::countHits(const BkgCluster& cluster, unsigned& nactive, unsigned& nstereo, const ComboHitCollection& chcol) const
  {
    nactive = nstereo = 0;
    for (const auto& chit : cluster.hits())
    {
      const ComboHit& ch = chcol[chit];
      nactive += ch.nStrawHits();
      if (ch.flag().hasAllProperties(stereo_)) nstereo += ch.nStrawHits();
    }
  }


  //----------------------------------------------------------------------------------
  int FlagBkgHits::findClusterIdx(BkgClusterCollection& bkgccol, unsigned ich) const
  {
    for (size_t icl=0;icl < bkgccol.size(); ++icl)
      if (std::find(bkgccol[icl].hits().begin(),bkgccol[icl].hits().end(),ich) != bkgccol[icl].hits().end()) return icl;
    return -1;
  }

}

DEFINE_ART_MODULE(mu2e::FlagBkgHits);
