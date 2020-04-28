//
// Two Niveau Threshold (TNT) algorithm, optimized version of two level threshold clustering originally developed from Dave Brown.
// Include fast preFilterting algorithm as well
//
//  Bertrand Echenard (2017) CIT
//
#ifndef TNTClusterer_HH
#define TNTClusterer_HH

#include "fhiclcpp/types/Atom.h"
#include "RecoDataProducts/inc/StrawDigiCollection.hh"
#include "TrkReco/inc/BkgClusterer.hh"
#include "fhiclcpp/types/Sequence.h"

#include <string>



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
              using Name = fhicl::Name;
              using Comment = fhicl::Comment;
              fhicl::Atom<float>            hitDistance {     Name("HitDistance"),      Comment("Minimum cluster hit distance")  };
              fhicl::Atom<float>            seedDistance{     Name("SeedDistance"),     Comment("Minimum distance for cluster seed")  };
              fhicl::Atom<float>            clusterDiameter{  Name("ClusterDiameter"),  Comment("Average cluster diameter")  };
              fhicl::Atom<float>            clusterTime{      Name("ClusterTime"),      Comment("Average cluster time spread")  };
              fhicl::Atom<float>            deltaTimeBinMin{  Name("DeltaTimeBinMin"),  Comment("Delta time for cluster lookup")  };
              fhicl::Atom<float>            maxHitTimeDiff{   Name("MaxHitTimeDiff"),   Comment("Maximum hit cluster tme difference")  };
              fhicl::Atom<float>            maxSumDistance{   Name("MaxSumDistance"),   Comment("Maximum sum pf hit-cluster distance for convergence") };        
              fhicl::Atom<float>            minHitError{      Name("MinHitError"),      Comment("Min value of hit error")  };
              fhicl::Atom<float>            maxDistance{      Name("MaxDistance"),      Comment("Max hit-cluster distance")  };
              fhicl::Atom<float>            timeRMS{          Name("TimeRMS"),          Comment("Cluster time RMS")  };
              fhicl::Atom<unsigned>         maxCluIterations{ Name("MaxCluIterations"), Comment("Maximum number of cluster algo iterations") };
              fhicl::Atom<bool>             medianCentroid {  Name("MedianCentroid"),   Comment("Use median to calculate cluster centroid") };
              fhicl::Atom<bool>             preFilter{        Name("preFilter"),        Comment("Fast preFiltering algorithm") };
              fhicl::Atom<float>            pfTimeBin{        Name("pfTimeBin"),        Comment("Time bin size for preFiltering algorithm") };
              fhicl::Atom<float>            pfPhiBin{         Name("pfPhiBin"),         Comment("Phi bin size for preFiltering algorithm") };
              fhicl::Atom<unsigned>         pfMinHit{         Name("pfMinHit"),         Comment("Minimum number of hits inside bin for preFitlering algorithm") };
              fhicl::Atom<unsigned>         pfMinSumHit{      Name("pfMinSumHit"),      Comment("Minimum number of hits for sum of bins for preFitlering algorithm") };
              fhicl::Atom<bool>             comboInit{        Name("ComboInit"),        Comment("Start with combo hits") };
              fhicl::Sequence<std::string>  bkgmsk{           Name("BackgroundMask"),   Comment("Bkg hit selection mask") };
              fhicl::Sequence<std::string>  sigmsk{           Name("SignalMask"),       Comment("Signal hit selection mask") };
              fhicl::Atom<bool>             testflag{         Name("TestFlag"),         Comment("Test hit flags") };
              fhicl::Atom<int>              diag{             Name("Diag"),             Comment("Diagnosis level"),0 };   
          };
          
          
          explicit TNTClusterer(const Config& config);
          virtual ~TNTClusterer() {};

          void          init();
          virtual void  findClusters(BkgClusterCollection& preFilterClusters, BkgClusterCollection& postFilterClusters, 
                                     const ComboHitCollection& shcol, float mbtime, int iev);
          virtual float distance(const BkgCluster& cluster, const ComboHit& hit) const; 


      private:                   
          static const int numBuckets = 256; //number of buckets to store the clusters vs time - optimized for speed
          using arrayVecBkg = std::array<std::vector<int>,numBuckets>;

          void     initClu(const ComboHitCollection& chcol, std::vector<BkgCluster>& clusters, std::vector<BkgHit>& hinfo, const std::vector<unsigned>& hitSel ); 
          void     preFilter(BkgClusterCollection& clusters, const ComboHitCollection& chcol, std::vector<unsigned>& hitSel, const float mbtime);

          void     clusterAlgo(const ComboHitCollection& chcol, std::vector<BkgCluster>& clusters, 
                               std::vector<BkgHit>& hinfo, float tbin);
          unsigned formClusters(const ComboHitCollection& chcol, std::vector<BkgCluster>& clusters, float tbin, 
                                std::vector<BkgHit>& hinfo, arrayVecBkg& hitIndex);
          void     mergeClusters(std::vector<BkgCluster>& clusters, const ComboHitCollection& chcol, 
                                 std::vector<BkgHit>& hinfo, float dt, float dd2);
          void     mergeTwoClu(BkgCluster& clu1, BkgCluster& clu2 );
          void     updateCluster(BkgCluster& cluster, const ComboHitCollection& chcol, std::vector<BkgHit>& hinfo);

          void     dump(const std::vector<BkgCluster>& clusters, std::vector<BkgHit>& hinfo);

          std::vector<int> hitDtIdx_;
          float            dhit_;      
          float            dseed_;      
          float            dd_;         
          float            dd2_; 
          float            dt_;
          float            maxwt_; 
          float            md2_;
          float            trms2inv_; 
          float            tbinMin_;        
          float            maxHitdt_;   
          float            maxDistSum_; 
          unsigned         maxNiter_;    
          bool             useMedian_;  
          bool             preFilter_;  
          float            pfTimeBin_;
          float            pfPhiBin_;
          unsigned         pfMinHit_;
          unsigned         pfMinSumHit_;
          bool             comboInit_;  
          StrawHitFlag     bkgmask_;    
          StrawHitFlag     sigmask_;    
	  bool	           testflag_;   
          int              diag_;
	  int              ditime_;
   };
}
#endif
