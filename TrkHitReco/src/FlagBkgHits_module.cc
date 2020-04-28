#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "ConditionsService/inc/AcceleratorParams.hh"

#include "MCDataProducts/inc/StrawDigiMCCollection.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/BkgCluster.hh"
#include "RecoDataProducts/inc/BkgClusterHit.hh"
#include "RecoDataProducts/inc/BkgQual.hh"

#include "Mu2eUtilities/inc/MVATools.hh"
#include "TrkReco/inc/TNTClusterer.hh"
#include "TrkReco/inc/ScanClusterer.hh"

#include <string>
#include <vector>


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
             fhicl::Atom<bool>                     filterOutput{         Name("FilterOutput"),         Comment("Produce filtered ComboHit collection")  };
             fhicl::Atom<bool>                     flagComboHits{        Name("FlagComboHits"),        Comment("Produce filtered flag comboHit collection") };
             fhicl::Atom<bool>                     flagStrawHits {       Name("FlagStrawHits"),        Comment("Produce filtered flag strawHit collection") };
             fhicl::Sequence<std::string>          backgroundMask{       Name("BackgroundMask"),       Comment("Bkg hit selection mask") };
             fhicl::Sequence<std::string>          stereoSelection{      Name("StereoSelection"),      Comment("Stereo hit selection mask") };
             fhicl::Atom<float>                    bkgMVAcut{            Name("BkgMVACut"),            Comment("Bkg MVA cut") };
             fhicl::Table<MVATools::Config>        bkgMVA{               Name("BkgMVA"),               Comment("MVA Configuration") };
             fhicl::Atom<bool>                     saveBkgClusters{      Name("SaveBkgClusters"),      Comment("Save bkg clusters"),false };
             fhicl::Atom<int>                      debugLevel{           Name("DebugLevel"),           Comment("Debug"),0 };
             fhicl::Atom<int>                      printFrequency{       Name("PrintFrequency"),       Comment("Print frequency"),100 };
             fhicl::Table<TNTClusterer::Config>    TNTClustering{        Name("TNTClustering"),        Comment("TNT Clusterer config") };
             fhicl::Table<ScanClusterer::Config>   ScanClustering{       Name("ScanClustering"),       Comment("Scan Clusterer config") };
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
         bool                                        filter_, flagch_, flagsh_;
         bool                                        savebkg_;
         StrawHitFlag                                bkgmsk_, stereo_;
         BkgClusterer*                               clusterer_;
         float                                       cperr2_;
         float                                       bkgMVAcut_;
         MVATools                                    bkgMVA_;
         int const                                   debug_;
         int const                                   printfreq_;
         int                                         iev_;

         void classifyCluster(BkgClusterCollection& bkgccolFast, BkgClusterCollection& bkgccol, BkgQualCollection& bkgqcol, 
                              StrawHitFlagCollection& chfcol, const ComboHitCollection& chcol) const;
         void fillBkgQual(    const BkgCluster& cluster, BkgQual& cqual, const ComboHitCollection& chcol) const;
         void fillMVA(        BkgQual& cqual) const;
         void countHits(      const BkgCluster& cluster, unsigned& nactive, unsigned& nstereo, const ComboHitCollection& chcol) const;
         void countPlanes(    const BkgCluster& cluster, BkgQual& cqual, const ComboHitCollection& chcol) const;
         int  findClusterIdx( BkgClusterCollection& bkgccol, unsigned ich) const;
  };


  FlagBkgHits::FlagBkgHits(const art::EDProducer::Table<Config>& config) :
     art::EDProducer{config},
     chtoken_{     consumes<ComboHitCollection>(config().comboHitCollection()) },
     shtoken_{     consumes<StrawHitCollection>(config().strawHitCollection()) },
     minnhits_(    config().minActiveHits() ),
     minnp_(       config().minNPlanes()),
     filter_(      config().filterOutput()),
     flagch_(      config().flagComboHits()),
     flagsh_(      config().flagStrawHits()),
     savebkg_(     config().saveBkgClusters()),
     bkgmsk_(      config().backgroundMask()),
     stereo_(      config().stereoSelection()),
     bkgMVAcut_(   config().bkgMVAcut()),
     bkgMVA_(      config().bkgMVA()),
     debug_(       config().debugLevel()),
     printfreq_(   config().printFrequency()),
     iev_(0)
  {
      // Must call consumesMany because fillStrawHitIndices calls getManyByType.
      consumesMany<ComboHitCollection>();

      if (flagch_) produces<StrawHitFlagCollection>("ComboHits");
      if (flagsh_) produces<StrawHitFlagCollection>("StrawHits");
      if (filter_) produces<ComboHitCollection>();
      if (savebkg_)
      {
          produces<BkgClusterHitCollection>();
          produces<BkgClusterCollection>();
          produces<BkgQualCollection>();
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
      bkgMVA_.initMVA();
  }


  //------------------------------------------------------------------------------------------
  void FlagBkgHits::produce(art::Event& event )
  {
    
     //unsigned iev=event.event();
      if (debug_ > 0 && iev_%printfreq_==0) std::cout<<"FlagBkgHits: event="<<iev_<<std::endl;

      ConditionsHandle<AcceleratorParams> accPar("ignored");
      float mbtime = accPar->deBuncherPeriod;

      auto chH = event.getValidHandle(chtoken_);
      const ComboHitCollection& chcol = *chH.product();
      unsigned nch = chcol.size();


      // the primary output is either a deep copy of selected inputs or a flag collection on those
      // intermediate results: keep these on the heap unless requested for diagnostics later
      BkgClusterCollection bkgccol,bkgccolFast;
      BkgQualCollection bkgqcol;
      BkgClusterHitCollection bkghitcol;
      bkgccol.reserve(nch/2);
      if (savebkg_) bkgqcol.reserve(bkgccol.size());
      if (savebkg_) bkghitcol.reserve(nch);


      // find clusters, sort is needed for recovery algorithm. bkgccolFast has hits that are autmoatically marked as bkg.
      clusterer_->findClusters(bkgccolFast,bkgccol,chcol, mbtime, iev_);
      sort(bkgccol.begin(),bkgccol.end(),[](const BkgCluster& c1,const BkgCluster& c2) {return c1.time() < c2.time();});


      // classify clusters                  
      StrawHitFlagCollection chfcol(nch);
      classifyCluster(bkgccolFast, bkgccol, bkgqcol, chfcol,chcol);


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
          event.put(std::make_unique<BkgQualCollection>(bkgqcol));
      }

      ++iev_;
      return;
  }
 
 
  //------------------------------------------------------------------------------------------
  void FlagBkgHits::classifyCluster(BkgClusterCollection& bkgccolFast, BkgClusterCollection& bkgccol, BkgQualCollection& bkgqcol, 
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
           BkgQual cqual;
           fillBkgQual(cluster, cqual, chcol);
           fillMVA(cqual);

           StrawHitFlag flag(StrawHitFlag::bkgclust);
           if (cqual.MVAOutput() > bkgMVAcut_)
           {
	      flag.merge(StrawHitFlag::bkg);
	      if (savebkg_) cluster._flag.merge(BkgClusterFlag::bkg);
           } 

           for (const auto& chit : cluster.hits()) chfcol[chit] = flag;
           if (savebkg_) bkgqcol.push_back(std::move(cqual));
      }
  }   
  
  
  //--------------------------------------------------------------------------------------------------------------
  void FlagBkgHits::fillBkgQual(const BkgCluster& cluster, BkgQual& cqual, const ComboHitCollection& chcol) const
  {
     unsigned nactive, nstereo;
     countHits(cluster, nactive, nstereo,chcol);
     cqual.setMVAStatus(MVAStatus::unset);

     if (nactive < minnhits_) return;
     countPlanes(cluster,cqual,chcol);
          
     cqual[BkgQual::hrho]  = 0;
     cqual[BkgQual::shrho] = 0;
     cqual[BkgQual::nhits] = nactive;
     cqual[BkgQual::sfrac] = static_cast<float>(nstereo)/nactive;
     cqual[BkgQual::crho] = sqrtf(cluster.pos().perp2());

     if (cqual[BkgQual::np] >= minnp_)
     {
         std::vector<float> hz;        
         for (const auto& chit : cluster.hits()) hz.push_back(chcol[chit].pos().z());

         // find the min, max and largest gap from the sorted Z positions
         std::sort(hz.begin(),hz.end());
         float zgap = 0.0;
         for (unsigned iz=1;iz<hz.size();++iz) zgap=std::max(zgap,hz[iz]-hz[iz-1]);
         cqual[BkgQual::zmin] = hz.front();
         cqual[BkgQual::zmax] = hz.back();
         cqual[BkgQual::zgap] = zgap;

         cqual.setMVAStatus(MVAStatus::filled);
     } 
     else 
     {
         cqual[BkgQual::zmin]  = -1.0;
         cqual[BkgQual::zmax]  = -1.0;
         cqual[BkgQual::zgap]  = -1.0;
     }      
  }
  
  
  //----------------------------------------------
  void FlagBkgHits::fillMVA(BkgQual& cqual) const
  {
       if (cqual.status() == MVAStatus::unset) return;

       std::vector<float> mvavars(7,0.0);
       mvavars[0] = cqual.varValue(BkgQual::crho);
       mvavars[1] = cqual.varValue(BkgQual::zmin);
       mvavars[2] = cqual.varValue(BkgQual::zmax);
       mvavars[3] = cqual.varValue(BkgQual::zgap);
       mvavars[4] = cqual.varValue(BkgQual::np);
       mvavars[5] = cqual.varValue(BkgQual::npfrac);
       mvavars[6] = cqual.varValue(BkgQual::nhits);
       float mvaout = bkgMVA_.evalMVA(mvavars);

       cqual.setMVAValue(mvaout);
       cqual.setMVAStatus(MVAStatus::calculated);
   }



   //-------------------------------------------------------------------------------------------------------------
   void FlagBkgHits::countPlanes(const BkgCluster& cluster, BkgQual& cqual, const ComboHitCollection& chcol) const 
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

        cqual[BkgQual::np]     = np;
        cqual[BkgQual::npexp]  = npexp;
        cqual[BkgQual::npfrac] = static_cast<float>(np)/static_cast<float>(npexp);
        cqual[BkgQual::nphits] = static_cast<float>(nphits)/static_cast<float>(np);
   }


   //----------------------------------------------------------------------------------------------------------------------------------
   void FlagBkgHits::countHits(const BkgCluster& cluster, unsigned& nactive, unsigned& nstereo, const ComboHitCollection& chcol) const 
   {
        nactive = nstereo = 0;
        for (const auto& chit : cluster.hits())
        {
           const ComboHit& ch = chcol[chit];
           nactive += ch.nStrawHits();
           if (ch.flag().hasAnyProperty(stereo_)) nstereo += ch.nStrawHits();
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
