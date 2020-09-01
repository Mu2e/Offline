//
// Fast clustering in Time-phi plane, followed by single pass clustering similar to TNTClusterer.
//
//  Bertrand Echenard (2020) CIT
//
#ifndef ScanClusterer_HH
#define ScanClusterer_HH

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "RecoDataProducts/inc/StrawDigiCollection.hh"
#include "TrkReco/inc/BkgClusterer.hh"


namespace mu2e {
 
   class ScanClusterer : public BkgClusterer
   {

      public:

         struct Config 
         {
              using Name    = fhicl::Name;
              using Comment = fhicl::Comment;
              fhicl::Atom<float>            tbin{             Name("Tbin"),             Comment("Time histo bin width") };
              fhicl::Atom<float>            pbin{             Name("Pbin"),             Comment("Phi histo bin width") };
              fhicl::Atom<float>            rbin{             Name("Rbin"),             Comment("Radius histo bin width") };
              fhicl::Atom<float>            rmin{             Name("Rmin"),             Comment("Min Radius histo bin width") };
              fhicl::Atom<float>            rmax{             Name("Rmax"),             Comment("Max Radius histo bin width") };
              fhicl::Atom<unsigned>         minPeakHit{       Name("MinPeakHit"),       Comment("Minimum number of hits in bin to be included") };
              fhicl::Atom<unsigned>         minSeedHit{       Name("MinSeedHit"),       Comment("Minimum number of hits in bin to start cluster") };
              fhicl::Atom<float>            minRadHit{        Name("MinRadHit"),        Comment("Radius difference to include hits") };           
              fhicl::Atom<unsigned>         filterAlgo{       Name("FilterAlgo"),       Comment("Filter Algorithm") };
              fhicl::Sequence<std::string>  bkgmsk{           Name("BackgroundMask"),   Comment("Bkg hit selection mask") };
              fhicl::Sequence<std::string>  sigmsk{           Name("SignalMask"),       Comment("Signal hit selection mask") };
              fhicl::Atom<bool>             testflag{         Name("TestFlag"),         Comment("Test hit flags") };
              fhicl::Atom<int>              diag{             Name("Diag"),             Comment("Diagnosis level"),0 };   
          };


          explicit ScanClusterer(const Config& config);
          virtual ~ScanClusterer() {};

          void          init();
          virtual void  findClusters(BkgClusterCollection& preFilterClusters, BkgClusterCollection& postFilterClusters, 
                                     const ComboHitCollection& shcol, float mbtime, int iev);

          virtual float distance(const BkgCluster& cluster, const ComboHit& hit) const {return 0.0;} 

      private:         
          
          void fastFilter1(BkgClusterCollection& clusters, const ComboHitCollection& chcol, const float mbtime);
          void fastFilter2(BkgClusterCollection& clusters, const ComboHitCollection& chcol, const float mbtime);
                

          float            tbin_; 
          float            pbin_;
          float            rbin_;
          float            rmin_;
          float            rmax_;
          unsigned         minPeakHit_;
          unsigned         minSeedHit_;
          float            minRadHit_;
          unsigned         filterAlgo_;
          StrawHitFlag     bkgmask_;    
          StrawHitFlag     sigmask_;    
	  bool	           testflag_;   
          int              diag_;
   };
}
#endif


