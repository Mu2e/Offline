// Clustering using Chi2 distance calculation method - Mete Yucel 2024
// Based on TNTClusterer

#ifndef Chi2Clusterer_HH
#define Chi2Clusterer_HH

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

#include "Offline/GeneralUtilities/inc/CombineTwoDPoints.hh"
#include "Offline/RecoDataProducts/inc/StrawDigi.hh"
#include "Offline/TrkHitReco/inc/BkgClusterer.hh"
#include "Offline/TrkHitReco/inc/TrainBkgDiag.hxx"
#include "Offline/TrkHitReco/inc/TrainBkgDiagStationChi2SLine.hxx"

#include <string>


//Inference class
namespace TMVA_SOFIE_TrainBkgDiagStationChi2SLine {
  class Session;
}


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
        fhicl::Atom<unsigned>         minActiveHits{    Name("MinActiveHits"),    Comment("Minumim number of active hits in a cluster") };
        fhicl::Atom<unsigned>         minNPlanes{       Name("MinNPlanes"),       Comment("Minumim number of planes in a cluster") };
        fhicl::Atom<std::string>      kerasWeights{     Name("KerasWeights"),     Comment("Weights for keras model") };
        fhicl::Atom<bool>             useSLine{         Name("UseSLine"),         Comment("Use SLine info") };
        fhicl::Atom<int>              diag{             Name("Diag"),             Comment("Diagnosis level"),0 };
      };

      Chi2Clusterer(const std::optional<Config> config);
      virtual ~Chi2Clusterer() {};

      void init();
      virtual void  findClusters   (BkgClusterCollection& clusters, const ComboHitCollection& shcol);
      virtual void  classifyCluster(BkgCluster& cluster,            const ComboHitCollection& chcol);
      virtual float distance       (const BkgCluster& cluster,      const ComboHit& hit) const;


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
      unsigned                minnhits_;
      unsigned                minnp_;
      std::string             kerasW_;
      bool                    useSLine_;
      BkgCluster::distMethod  distMethodFlag_;
      int                     diag_;

      std::shared_ptr<TMVA_SOFIE_TrainBkgDiag::Session>                 sofiePtr1_;
      std::shared_ptr<TMVA_SOFIE_TrainBkgDiagStationChi2SLine::Session> sofiePtr2_;
  };
}
#endif
