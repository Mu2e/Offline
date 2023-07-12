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
#include "Offline/TrkHitReco/inc/TrainBkgDiag.hxx"

#include <string>
#include <vector>

//Inference class
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
        fhicl::Atom<std::string>              kerasWeights{         Name("KerasWeights"),         Comment("Weights for keras model") };
        fhicl::Atom<float>                    kerasQuality{         Name("KerasQuality"),         Comment("Keras quality cut") };
      };

      enum clusterer {TwoNiveauThreshold=1};
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
      std::string                                 kerasW_;
      float                                       kerasQ_;
      int                                         iev_;
      std::shared_ptr<TMVA_SOFIE_TrainBkgDiag::Session> sofiePtr;

      void classifyCluster(BkgClusterCollection& bkgccol, StrawHitFlagCollection& chfcol, const ComboHitCollection& chcol) const;
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
    kerasW_(      config().kerasWeights()),
    kerasQ_(      config().kerasQuality()),
    iev_(0)
    {
      // Must call consumesMany because fillStrawHitIndices calls getManyByType.
      consumesMany<ComboHitCollection>();
      ConfigFileLookupPolicy configFile;
      auto kerasWgtsFile = configFile(kerasW_);
      sofiePtr = std::make_shared<TMVA_SOFIE_TrainBkgDiag::Session>(kerasWgtsFile);

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
    BkgClusterCollection bkgccol;
    BkgClusterHitCollection bkghitcol;
    bkgccol.reserve(nch/2);
    if (savebkg_) bkghitcol.reserve(nch);


    // find clusters, sort is needed for recovery algorithm. bkgccolFast has hits that are autmoatically marked as bkg.
    clusterer_->findClusters(bkgccol,chcol, iev_);
    std::sort(bkgccol.begin(),bkgccol.end(),[](const BkgCluster& c1,const BkgCluster& c2) {return c1.time() < c2.time();});

    // classify clusters
    StrawHitFlagCollection chfcol(nch);
    classifyCluster(bkgccol, chfcol, chcol);

    //produce BkgClusterHit info collection
    if (savebkg_) {
      for (size_t ich=0;ich < chcol.size(); ++ich) {
        const ComboHit& ch = chcol[ich];
        int icl = findClusterIdx(bkgccol,ich);
        if (icl > -1) bkghitcol.emplace_back(BkgClusterHit(clusterer_->distance(bkgccol[icl],ch),ch.flag()));
        else          bkghitcol.emplace_back(BkgClusterHit(999.0,ch.flag()));
      }
    }

    //produce filtered ComboHit collection
    if (filter_) {
      auto chcolFilter = std::make_unique<ComboHitCollection>();
      chcolFilter->reserve(nch);
      // same parent as the original collection
      chcolFilter->setParent(chcol.parent());
      for (size_t ich=0;ich < nch; ++ich) {
        StrawHitFlag const& flag = chfcol[ich];
        if (!flag.hasAnyProperty(bkgmsk_)) {
          chcolFilter->push_back(chcol[ich]);
          chcolFilter->back()._flag.merge(flag);
        }
      }
      event.put(std::move(chcolFilter));
    }

    //produce StrawHit flags
    // Note: fillStrawHitIndices is quite slow and should be improved
    if (flagsh_) {
      auto shH = event.getValidHandle(shtoken_);
      const StrawHitCollection* shcol  = shH.product();

      unsigned nsh = shcol->size();
      auto shfcol  = std::make_unique<StrawHitFlagCollection>(nsh);
      std::vector<std::vector<StrawHitIndex> > shids;
      chcol.fillStrawHitIndices(event,shids);
      for (size_t ich = 0;ich < nch;++ich) {
        StrawHitFlag flag = chfcol[ich];
        flag.merge(chcol[ich].flag());
        for(auto ish : shids[ich]) (*shfcol)[ish] = flag;
      }

      event.put(std::move(shfcol),"StrawHits");
    }

    //produce ComboHit flags
    if (flagch_) {
      for(size_t ich=0;ich < nch; ++ich) chfcol[ich].merge(chcol[ich].flag());
      event.put(std::make_unique<StrawHitFlagCollection>(std::move(chfcol)),"ComboHits");
    }

    //produce background collection
    if (savebkg_) {
      event.put(std::make_unique<BkgClusterHitCollection>(bkghitcol));
      event.put(std::make_unique<BkgClusterCollection>(bkgccol));
    }

    ++iev_;
    return;
  }


  //------------------------------------------------------------------------------------------
  void FlagBkgHits::classifyCluster(BkgClusterCollection& bkgccol, StrawHitFlagCollection& chfcol, const ComboHitCollection& chcol) const
  {
    for (auto& cluster : bkgccol) {
      // count hits and planes
      std::array<int,StrawId::_nplanes> hitplanes{0};
      for (const auto& chit : cluster.hits()) {
        const ComboHit& ch = chcol[chit];
        hitplanes[ch.strawId().plane()] += ch.nStrawHits();
      }
      unsigned npexp(0),np(0),nhits(0);
      int ipmin(0),ipmax(StrawId::_nplanes-1);
      while (hitplanes[ipmin]==0 && ipmin<StrawId::_nplanes) ++ipmin;
      while (hitplanes[ipmax]==0 and ipmax>0)                --ipmax;
      int fp(ipmin),lp(ipmin-1),pgap(0);
      for (int ip = ipmin; ip <= ipmax; ++ip) {
        npexp++; // should use TTracker to see if plane is physically present FIXME!
        if (hitplanes[ip]> 0){
          ++np;
          if(lp > 0 && ip - lp -1 > pgap)pgap = ip - lp -1;
          if(ip > lp)lp = ip;
          if(ip < fp)fp = ip;
          lp = ip;
        }
        nhits += hitplanes[ip];
      }
      if(nhits >= minnhits_ && np >= minnp_){
        // find averages
        double sumEdep(0.);
        double sqrSumDeltaTime(0.);
        double sqrSumDeltaX(0.);
        double sqrSumDeltaY(0.);
        unsigned nchits = cluster.hits().size();
        for (const auto& chit : cluster.hits()) {
          sumEdep +=  chcol[chit].energyDep()/chcol[chit].nStrawHits();
          sqrSumDeltaX += std::pow(chcol[chit].pos().x() - cluster.pos().x(),2);
          sqrSumDeltaY += std::pow(chcol[chit].pos().y() - cluster.pos().y(),2);
          sqrSumDeltaTime += std::pow(chcol[chit].time() - cluster.time(),2);
        }
        // fill mva input variables
        std::array<float,9> kerasvars;
        kerasvars[0] = cluster.pos().Rho(); // cluster rho, cyl coor
        kerasvars[1] = fp;// first plane hit
        kerasvars[2] = lp;// last plane hit
        kerasvars[3] = pgap;// largest plane gap without hits between planes with hits
        kerasvars[4] = np;// # of planes hit
        kerasvars[5] =  static_cast<float>(np)/static_cast<float>(lp - fp);// fraction of planes hit between first and last plane
        kerasvars[6] = nhits;// sum of straw hits
        kerasvars[7] = std::sqrt((sqrSumDeltaX+sqrSumDeltaY)/nchits);  // RMS of cluster rho
        kerasvars[8] = std::sqrt(sqrSumDeltaTime/nchits);// RMS of cluster time

        auto kerasout = sofiePtr->infer(kerasvars.data());
        cluster.setKerasQ(kerasout[0]);
        if(debug_>0)std::cout << "kerasout = " << kerasout[0] << std::endl;

        StrawHitFlag flag(StrawHitFlag::bkgclust);
        if (cluster.getKerasQ()> kerasQ_) {
          flag.merge(StrawHitFlag(StrawHitFlag::bkg));
          cluster._flag.merge(BkgClusterFlag::bkg);
        }
        for (const auto& chit : cluster.hits()) chfcol[chit].merge(flag);
      } else
        cluster.setKerasQ(-1.0);
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

DEFINE_ART_MODULE(mu2e::FlagBkgHits)
