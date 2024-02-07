// Clustering using Chi2 distance calculation method - Mete Yucel 2024
// Based on TNTClusterer

#ifndef Chi2Clusterer_HH
#define Chi2Clusterer_HH

#include "fhiclcpp/types/Atom.h"
#include "Offline/RecoDataProducts/inc/StrawDigi.hh"
#include "Offline/TrkHitReco/inc/BkgClusterer.hh"
#include "Offline/GeneralUtilities/inc/CombineTwoDPoints.hh"
#include "fhiclcpp/types/Sequence.h"

#include <string>



namespace mu2e {

  struct Chi2BkgHit
  {
    Chi2BkgHit(unsigned chidx): chidx_(chidx),clusterIdx_(-1) {};

    unsigned     chidx_;
    int          clusterIdx_;

  };


  class Chi2Clusterer : public BkgClusterer
  {
    public:

      struct Config
      {
        using Name    = fhicl::Name;
        using Comment = fhicl::Comment;
        fhicl::Atom<float>            hitDistance{      Name("HitDistance"),      Comment("Minimum cluster hit distance")  };
        fhicl::Atom<float>            seedDistance{     Name("SeedDistance"),     Comment("Minimum distance for cluster seed")  };
        fhicl::Atom<unsigned>         minClusterHits{   Name("MinClusterHits"),   Comment("Cut for minimum cluster hit size"),1 };
        fhicl::Atom<float>            maxHitTimeDiff{   Name("MaxHitTimeDiff"),   Comment("Maximum hit cluster time difference")  };
        fhicl::Atom<unsigned>         maxCluIterations{ Name("MaxCluIterations"), Comment("Maximum number of cluster algo iterations") };
        fhicl::Sequence<std::string>  bkgmsk{           Name("BackgroundMask"),   Comment("Bkg hit selection mask") };
        fhicl::Sequence<std::string>  sigmsk{           Name("SignalMask"),       Comment("Signal hit selection mask") };
        fhicl::Atom<bool>             testflag{         Name("TestFlag"),         Comment("Test hit flags") };
        fhicl::Atom<int>              diag{             Name("Diag"),             Comment("Diagnosis level"),0 };
      };

      Chi2Clusterer(const std::optional<Config> config);
      virtual ~Chi2Clusterer() {};

      void init();
      virtual void  findClusters(BkgClusterCollection& clusters, const ComboHitCollection& shcol, int iev);
      virtual float distance    (const BkgCluster& cluster,      const ComboHit& hit) const;


    private:
      size_t numBuckets_; //number of buckets to store the cluster ids vs time

      void     initClustering  (const ComboHitCollection& chcol, std::vector<Chi2BkgHit>& hinfo);
      void     doClustering    (const ComboHitCollection& chcol, std::vector<BkgCluster>& clusters, std::vector<Chi2BkgHit>& hinfo);
      unsigned formClusters    (const ComboHitCollection& chcol, std::vector<BkgCluster>& clusters, std::vector<Chi2BkgHit>& hinfo);
      void     dump            (const std::vector<BkgCluster>& clusters, const std::vector<Chi2BkgHit>& hinfo);

      float                   tbin_;
      float                   tmin_;
      float                   tmax_;
      float                   dhit_;
      float                   dseed_;
      float                   dt_;
      float                   chi2Cut_;
      unsigned                minClusterHits_;
      float                   maxHitdt_;
      unsigned                maxNiter_;
      StrawHitFlag            bkgmask_;
      StrawHitFlag            sigmask_;
      bool                    testflag_;
      int                     diag_;
      BkgCluster::distMethod  distMethodFlag_;
  };
}
#endif
