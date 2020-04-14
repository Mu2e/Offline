//
// Fast clustering in Time-phi plane, followed by single pass clustering similar to TNTClusterer.
//
//  Bertrand Echenard (2020) CIT
//
#ifndef ScanTPClusterer_HH
#define ScanTPClusterer_HH

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "RecoDataProducts/inc/StrawDigiCollection.hh"
#include "TrkReco/inc/BkgClusterer.hh"


namespace mu2e {


   class ScanTPClusterer : public BkgClusterer
   {

      public:

         struct Config 
         {
              using Name    = fhicl::Name;
              using Comment = fhicl::Comment;
              fhicl::Atom<float>            tmin{             Name("Tmin"),             Comment("Minimum hit time") };
              fhicl::Atom<float>            tmax{             Name("Tmax"),             Comment("Maximum hit time") };
              fhicl::Atom<float>            tbin{             Name("Tbin"),             Comment("Time histo bin width") };
              fhicl::Atom<float>            pbin{             Name("Pbin"),             Comment("Phi histo bin width") };
              fhicl::Atom<unsigned>         minPeakHit{       Name("MinPeakHit"),       Comment("Minimum number of hits in bin to be included") };
              fhicl::Atom<unsigned>         minSeedHit{       Name("MinSeedHit"),       Comment("Minimum number of hits in bin to start cluster") };
              fhicl::Atom<float>            minRadHit{        Name("MinRadHit"),        Comment("Radius difference to include hits") };           
              fhicl::Atom<float>            hitDistance{      Name("HitDistance"),      Comment("Minimum cluster hit distance") };
              fhicl::Atom<float>            seedDistance{     Name("SeedDistance"),     Comment("Minimum distance for cluster seed") };
              fhicl::Atom<float>            clusterDiameter{  Name("ClusterDiameter"),  Comment("Average cluster diameter") };
              fhicl::Atom<float>            clusterTime{      Name("ClusterTime"),      Comment("Average cluster time spread") };
              fhicl::Atom<float>            deltaTimeBinMin{  Name("DeltaTimeBinMin"),  Comment("Delta time for cluster lookup") };
              fhicl::Atom<float>            maxHitTimeDiff{   Name("MaxHitTimeDiff"),   Comment("Maximum hit cluster tme difference") };              
              fhicl::Atom<float>            minHitError{      Name("MinHitError"),      Comment("Min value of hit error") };
              fhicl::Atom<float>            maxDistance{      Name("MaxDistance"),      Comment("Max hit-cluster distance") };
              fhicl::Atom<float>            timeRMS{          Name("TimeRMS"),          Comment("Cluster time RMS") };
              fhicl::Atom<unsigned>         filterAlgo{       Name("FilterAlgo"),       Comment("Fast filter algorithm") };
              fhicl::Atom<bool>             comboInit{        Name("ComboInit"),        Comment("Start with combo hits") };
              fhicl::Sequence<std::string>  bkgmsk{           Name("BackgroundMask"),   Comment("Bkg hit selection mask") };
              fhicl::Sequence<std::string>  sigmsk{           Name("SignalMask"),       Comment("Signal hit selection mask") };
              fhicl::Atom<bool>             testflag{         Name("TestFlag"),         Comment("Test hit flags") };
              fhicl::Atom<int>              diag{             Name("Diag"),             Comment("Diagnosis level"),0 };   
          };


          explicit ScanTPClusterer(const Config& config);
          virtual ~ScanTPClusterer() {};

          void          init();
          virtual void  findClusters(BkgClusterCollection& preFilterClusters, BkgClusterCollection& postFilterClusters, 
                                     const ComboHitCollection& shcol, float mbtime, int iev);
          virtual float distance(const BkgCluster& cluster, const ComboHit& hit) const; 


      private:         
          
          static const int numBuckets = 200; //number of buckets to store the clusters vs time - seems fine
          void fastFilter1(BkgClusterCollection& clusters, const ComboHitCollection& chcol, std::vector<unsigned>& hitSel);
          void fastFilter2(BkgClusterCollection& clusters, const ComboHitCollection& chcol, std::vector<unsigned>& hitSel);
          void fastCluster(BkgClusterCollection& clusters, const ComboHitCollection& chcol, const std::vector<unsigned>& hitSel, float mbtime);
       

          float            tmin_, tmax_, tbin_, pbin_;
          unsigned         minPeakHit_, minSeedHit_;
          float            minRadHit_;
          float            dhit_;      
          float            dseed_;      
          float            dd_;         
          float            dt_;
          float            tbinMinCluster_;        
          float            maxHitdt_;   
          unsigned         filterAlgo_;
          bool             comboInit_; 
          bool             addOnePassClu_; 
          StrawHitFlag     bkgmask_;    
          StrawHitFlag     sigmask_;    
	  bool	           testflag_;   
          int              diag_;
          float            trms2inv_; 
	  int              ditime_;
          float            dd2_; 
          float            maxwt_; 
          float            md2_;
   };
}
#endif
