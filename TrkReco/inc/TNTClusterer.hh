//
// Two Niveau Threshold (TNT) algorithm, optimized version of two level threshold clustering originally developed from Dave Brown.
//
//  Bertrand Echenard (2017) CIT
//
#ifndef TNTClusterer_HH
#define TNTClusterer_HH

#include "RecoDataProducts/inc/StrawDigiCollection.hh"
#include "TrkReco/inc/BkgClusterer.hh"
#include "TrkReco/inc/ClustererFclConfig.hh"


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
          explicit TNTClusterer(const TNTClustererConfig& config);
          virtual ~TNTClusterer() {};

          void          init();
          virtual void  findClusters(BkgClusterCollection& clusters, const ComboHitCollection& shcol, float mbtime, int iev);
          virtual float distance(const BkgCluster& cluster, const ComboHit& hit) const; 


      private:                   
          static const int numBuckets = 256; //number of buckets to store the clusters vs time - optimized for speed
          using arrayVecBkg = std::array<std::vector<int>,numBuckets>;

          void     initClu(const ComboHitCollection& chcol, std::vector<BkgCluster>& clusters, std::vector<BkgHit>& hinfo); 
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
          bool             comboInit_;  
          StrawHitFlag     bkgmask_;    
          StrawHitFlag     sigmask_;    
	  bool	           testflag_;   
          int              diag_;
	  int              ditime_;
   };
}
#endif
