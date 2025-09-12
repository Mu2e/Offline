//
// Two Niveau Threshold (TNT) algorithm, optimized version of two level threshold clustering originally developed from Dave Brown.
// Include fast preFilterting algorithm as well
//
//  Bertrand Echenard (2017) CIT
//
#ifndef TNTClusterer_HH
#define TNTClusterer_HH

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "Offline/RecoDataProducts/inc/StrawDigi.hh"
#include "Offline/TrkHitReco/inc/BkgClusterer.hh"
#include "Offline/TrkHitReco/inc/TrainBkgDiag.hxx"

#include <string>


//Inference class
namespace TMVA_SOFIE_TrainBkgDiag {
  class Session;
}


namespace mu2e {

  struct BkgHit
  {
    BkgHit(unsigned chidx): chidx_(chidx),distance_(1000.0),clusterIdx_(-1) {};

    unsigned     chidx_;
    float        distance_;
    int          clusterIdx_;
  };


  class TNTClusterer : public BkgClusterer
  {
    public:

      struct Config
      {
        using Name    = fhicl::Name;
        using Comment = fhicl::Comment;
        fhicl::Atom<float>            hitDistance{      Name("HitDistance"),      Comment("Minimum cluster hit distance")  };
        fhicl::Atom<float>            seedDistance{     Name("SeedDistance"),     Comment("Minimum distance for cluster seed")  };
        fhicl::Atom<float>            clusterDiameter{  Name("ClusterDiameter"),  Comment("Average cluster diameter")  };
        fhicl::Atom<float>            clusterTime{      Name("ClusterTime"),      Comment("Average cluster time spread")  };
        fhicl::Atom<int>              minClusterHits{   Name("MinClusterHits"),   Comment("Cut for minimum cluster hit size"),1 };
        fhicl::Atom<float>            maxHitTimeDiff{   Name("MaxHitTimeDiff"),   Comment("Maximum hit cluster tme difference")  };
        fhicl::Atom<float>            maxSumDistance{   Name("MaxSumDistance"),   Comment("Maximum sum pf hit-cluster distance for convergence") };
        fhicl::Atom<float>            minHitError{      Name("MinHitError"),      Comment("Min value of hit error")  };
        fhicl::Atom<float>            maxDistance{      Name("MaxDistance"),      Comment("Max hit-cluster distance")  };
        fhicl::Atom<float>            timeRMS{          Name("TimeRMS"),          Comment("Cluster time RMS")  };
        fhicl::Atom<unsigned>         maxCluIterations{ Name("MaxCluIterations"), Comment("Maximum number of cluster algo iterations") };
        fhicl::Atom<bool>             medianCentroid {  Name("MedianCentroid"),   Comment("Use median to calculate cluster centroid") };
        fhicl::Atom<bool>             comboInit{        Name("ComboInit"),        Comment("Start with combo hits") };
        fhicl::Sequence<std::string>  bkgmsk{           Name("BackgroundMask"),   Comment("Bkg hit selection mask") };
        fhicl::Sequence<std::string>  sigmsk{           Name("SignalMask"),       Comment("Signal hit selection mask") };
        fhicl::Atom<bool>             testflag{         Name("TestFlag"),         Comment("Test hit flags") };
        fhicl::Atom<unsigned>         minActiveHits{    Name("MinActiveHits"),    Comment("Minumim number of active hits in a cluster") };
        fhicl::Atom<unsigned>         minNPlanes{       Name("MinNPlanes"),       Comment("Minumim number of planes in a cluster") };
        fhicl::Atom<std::string>      kerasWeights{     Name("KerasWeights"),     Comment("Weights for keras model") };
        fhicl::Atom<int>              diag{             Name("Diag"),             Comment("Diagnosis level"),0 };
      };


      TNTClusterer(const std::optional<Config> config);
      virtual ~TNTClusterer() {};

      void          init        ();
      virtual void  findClusters   (BkgClusterCollection& clusters, const ComboHitCollection& shcol);
      virtual void  classifyCluster(BkgCluster& cluster,            const ComboHitCollection& chcol);
      virtual float distance       (const BkgCluster& cluster,      const ComboHit& hit) const;


    private:
      static constexpr int numBuckets_ =256; //number of buckets to store the cluster ids vs time

      void     initClustering  (const ComboHitCollection& chcol, std::vector<BkgHit>& hinfo);
      void     doClustering    (const ComboHitCollection& chcol, std::vector<BkgCluster>& clusters, std::vector<BkgHit>& hinfo);
      unsigned formClusters    (const ComboHitCollection& chcol, std::vector<BkgCluster>& clusters, std::vector<BkgHit>& hinfo);
      void     mergeClusters   (std::vector<BkgCluster>& clusters, const ComboHitCollection& chcol, std::vector<BkgHit>& hinfo,
                                float dt, float dd2);
      void     mergeTwoClusters(BkgCluster& clu1, BkgCluster& clu2);
      void     updateCluster   (BkgCluster& cluster, const ComboHitCollection& chcol, std::vector<BkgHit>& hinfo);
      void     dump            (const std::vector<BkgCluster>& clusters, const std::vector<BkgHit>& hinfo);

      float                   tbin_;
      float                   dhit_;
      float                   dseed_;
      float                   dd_;
      float                   dd2_;
      float                   dt_;
      int                     minClusterHits_;
      float                   maxwt_;
      float                   md2_;
      float                   trms2inv_;
      float                   maxHitdt_;
      float                   maxDistSum_;
      unsigned                maxNiter_;
      bool                    useMedian_;
      bool                    comboInit_;
      StrawHitFlag            bkgmask_;
      StrawHitFlag            sigmask_;
      bool                    testflag_;
      unsigned                minnhits_;
      unsigned                minnp_;
      std::string             kerasW_;
      int                     diag_;
      BkgCluster::distMethod  distMethodFlag_;

      std::shared_ptr<TMVA_SOFIE_TrainBkgDiag::Session> sofiePtr_;
  };
}
#endif
