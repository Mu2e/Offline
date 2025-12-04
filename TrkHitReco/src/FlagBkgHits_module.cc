#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/OptionalTable.h"

#include "Offline/ConditionsService/inc/ConditionsHandle.hh"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "Offline/DataProducts/inc/StrawIdMask.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/BkgCluster.hh"
#include "Offline/RecoDataProducts/inc/BkgClusterHit.hh"
#include "Offline/TrkHitReco/inc/TNTClusterer.hh"
#include "Offline/TrkHitReco/inc/Chi2Clusterer.hh"
#include "Offline/TrkHitReco/inc/DBSClusterer.hh"

#include <string>
#include <vector>


namespace mu2e
{

  class FlagBkgHits : public art::EDProducer
  {
    public:

      struct Config
      {
        using Name    = fhicl::Name;
        using Comment = fhicl::Comment;

        fhicl::Atom<art::InputTag>                    comboHitCollection{   Name("ComboHitCollection"),   Comment("ComboHit collection name") };
        fhicl::Atom<float>                            clusterPositionError{ Name("ClusterPositionError"), Comment("Cluster poisiton error") };
        fhicl::Atom<int>                              clusterAlgorithm{     Name("ClusterAlgorithm"),     Comment("Clusterer algorithm") };
        fhicl::Atom<bool>                             filterHits{           Name("FilterHits"),           Comment("Produce filtered ComboHit collection")  };
        fhicl::Sequence<std::string>                  backgroundMask{       Name("BackgroundMask"),       Comment("Bkg hit selection mask") };
        fhicl::Atom<bool>                             saveBkgClusters{      Name("SaveBkgClusters"),      Comment("Save bkg clusters") };
        fhicl::Atom<std::string>                      outputLevel{          Name("OutputLevel"),          Comment("Level of the output ComboHitCollection") };
        fhicl::OptionalTable<TNTClusterer::Config>    TNTClustering{        Name("TNTClustering"),        Comment("TNT Clusterer config") };
        fhicl::OptionalTable<Chi2Clusterer::Config>   Chi2Clustering{       Name("Chi2Clustering"),       Comment("Chi2 Clusterer config") };
        fhicl::OptionalTable<DBSClusterer::Config>    DBSClustering{        Name("DBSClustering"),        Comment("DBS Clusterer config") };
        fhicl::Atom<float>                            kerasQuality{         Name("KerasQuality"),         Comment("Keras quality cut") };
        fhicl::Atom<int>                              debugLevel{           Name("DebugLevel"),           Comment("Debug"),0 };
      };

      enum clusterer {TNT=1, Chi2=2, DBS=3};
      explicit FlagBkgHits(const art::EDProducer::Table<Config>& config);
      void beginJob() override;
      void produce(art::Event& event) override;

    private:
      const art::ProductToken<ComboHitCollection> chtoken_;
      bool                                        filter_;
      bool                                        savebkg_;
      StrawHitFlag                                bkgmsk_;
      StrawIdMask::Level                          level_;
      std::unique_ptr<BkgClusterer>               clusterer_;
      float                                       cperr2_;
      int const                                   debug_;
      float                                       kerasQ_;
      int                                         iev_;

      void classifyCluster(BkgClusterCollection& bkgccol, StrawHitFlagCollection& chfcol, const ComboHitCollection& chcol) const;
      int  findClusterIdx( BkgClusterCollection& bkgccol, unsigned ich) const;
  };



  FlagBkgHits::FlagBkgHits(const art::EDProducer::Table<Config>& config) :
    art::EDProducer{config},
    chtoken_{     consumes<ComboHitCollection>(config().comboHitCollection()) },
    filter_(      config().filterHits()),
    savebkg_(     config().saveBkgClusters()),
    bkgmsk_(      config().backgroundMask()),
    debug_(       config().debugLevel()),
    kerasQ_(      config().kerasQuality()),
    iev_(0)
    {
      ConfigFileLookupPolicy configFile;
      produces<ComboHitCollection>();

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
        case TNT:
          if (!config().TNTClustering())
          {
            throw cet::exception("RECO")<< "FlagBkgHits: TNTClusterer is not configured. Configure by adding\n"
                                        << "physics.producers.FlagBkgHits.TNTClustering : {@table::TNTClusterer}" << std::endl;
          }
          clusterer_ = std::make_unique<TNTClusterer>(config().TNTClustering());
          break;

        case Chi2:
          if (!config().Chi2Clustering())
          {
            throw cet::exception("RECO")<< "FlagBkgHits: Chi2Clusterer is not configured. Configure by adding\n"
                                        << "physics.producers.FlagBkgHits.Chi2Clustering : {@table::Chi2Clusterer}" << std::endl;
          }
          clusterer_ = std::make_unique<Chi2Clusterer>(config().Chi2Clustering());
          break;

       case DBS:
          if (!config().DBSClustering())
          {
            throw cet::exception("RECO") << "FlagBkgHits: DBSClusterer is not configured. Configure by adding\n"
                                        << "physics.producers.FlagBkgHits.DBSClustering : {@table::DBSClusterer}" << std::endl;
          }
          clusterer_ = std::make_unique<DBSClusterer>(config().DBSClustering());
          break;

        default:
          throw cet::exception("RECO")<< "Unknown clusterer" << ctype << std::endl;
      }

      StrawIdMask mask(config().outputLevel());
      level_ = mask.level();
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
    clusterer_->findClusters(bkgccol, chcol);
    std::sort(bkgccol.begin(),bkgccol.end(),[](const BkgCluster& c1,const BkgCluster& c2) {return c1.time() < c2.time();});

    // classify clusters
    StrawHitFlagCollection chfcol(nch);
    classifyCluster(bkgccol, chfcol, chcol);

    //produce BkgClusterHit info collection
    if (savebkg_) {
      for (size_t ich=0;ich < nch; ++ich) {
        const ComboHit& ch = chcol[ich];
        //StrawHitFlag const& flagbkg = chfcol[ich];
        int icl = findClusterIdx(bkgccol,ich);
        if (icl > -1) bkghitcol.emplace_back(BkgClusterHit(clusterer_->distance(bkgccol[icl],ch),ch.flag()));
        else          bkghitcol.emplace_back(BkgClusterHit(999.0,ch.flag()));
      }
    }

    // produce output ComboHit collection, either filtered or not

    auto chcol_out = std::make_unique<ComboHitCollection>();
    if(chfcol.size()>0){
      // same parent as the original collection
      if(level_ == chcol.level()){
        chcol_out->setSameParent(chcol);
        chcol_out->reserve(nch);
        for(size_t ich=0;ich < nch; ++ich) {
          StrawHitFlag const& flag = chfcol[ich];
          //std::cout<<"Before Flag = "<<flag<<std::endl;
          if (! filter_ || !flag.hasAnyProperty(bkgmsk_)) {
            // write out hits
            chcol_out->push_back(chcol[ich]);
            chcol_out->back()._flag.merge(flag);
          }
        }
      } else {
        // go down to the specified level
        auto pptr = chcol.parent(level_);
        auto const& chcol_p = *pptr;
        chcol_out->setSameParent(chcol_p);
        ComboHitCollection::SHIV shiv;
        chcol_out->reserve(chcol_p.size()*2);
        for(size_t ich=0;ich < nch; ++ich) {
          shiv.clear();
          StrawHitFlag const& flag = chfcol[ich];
          if (! filter_ || !flag.hasAnyProperty(bkgmsk_)) {
            // write out hits
            if(&chcol_p == chcol.fillStrawHitIndices(ich,shiv,level_)){
              for(auto ishi : shiv ){
                chcol_out->push_back(chcol_p[ishi]);
                chcol_out->back()._flag.merge(flag);
              }
            } else {
              throw cet::exception("RECO")<< "FlagBkgHits: inconsistent ComboHits" << std::endl;
            }
          }
        }
        if((! filter_) && chcol_out->size() != chcol_p.size())
          throw cet::exception("RECO")<< "FlagBkgHits: inconsistent ComboHit output" << std::endl;
      }
    }
    event.put(std::move(chcol_out));

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
      clusterer_->classifyCluster(cluster,chcol);
      StrawHitFlag flag(StrawHitFlag::bkgclust);
      if (cluster.getKerasQ()> kerasQ_) {
        flag.merge(StrawHitFlag(StrawHitFlag::bkg));
        cluster._flag.merge(BkgClusterFlag::bkg);
      }

      for (const auto& chit : cluster.hits())
        chfcol[chit].merge(flag);
    }
    /*for (auto& cluster : bkgccol) {
        std::cout << "Cluster has " << cluster.hitposition().size() << " hit positions\n";
        for (size_t i = 0; i < cluster.hitposition().size(); ++i) {
          auto const& p = cluster.hitposition()[i];
          std::cout << "  hit " << i
                  << " pos = (" << p.x()
                  << ", " << p.y()
                  << ", " << p.z()
                  << ")\n";
        }
        }*/
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
