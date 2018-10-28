//
// Two Niveau algorithm, optimized version of TLT 
//
//  Bertrand Echenard (2017) CIT
//
#ifndef TNTClusterer_HH
#define TNTClusterer_HH

#include "TrkReco/inc/BkgClusterer.hh"
#include "DataProducts/inc/XYZVec.hh"
#include "fhiclcpp/ParameterSet.h"
#include "TTree.h"
#include "TH2.h"
#include "TH1.h"

#include "MCDataProducts/inc/StrawDigiMCCollection.hh"
#include "RecoDataProducts/inc/StrawDigiCollection.hh"

namespace mu2e {


  class ClusterStraw {  
  
     public:

        ClusterStraw(BkgClusterHit& hit, const ComboHit& chit);
        virtual ~ClusterStraw() {};

        inline const  XYZVec& pos()                       const { return _pos; }
        inline float  time()                              const { return _time; }
        inline int    itime()                             const { return int(_time/10.0); }
        inline int    id()                                const { return _id; }
        inline bool   hasChanged()                        const { return _hasChanged;}
        inline void   flagChanged(bool val)                     { _hasChanged = val;}
        inline const  std::vector<BkgClusterHit*>& hits() const { return _hitsPtr; }
        inline        std::vector<BkgClusterHit*>& hits()       { return _hitsPtr; }
        static void   resetCounter()                            { idCounter_ = 0; }

        void updateCache(const ComboHitCollection& chcol, float _maxwt);
  
  
     private:
	int                         _id;
        XYZVec                      _pos;
        float                       _time;
        bool                        _hasChanged;
        std::vector<BkgClusterHit*> _hitsPtr; 

        static int idCounter_;  
  };



    
  class TNTClusterer : public BkgClusterer
  {
     
     public:

         explicit TNTClusterer(fhicl::ParameterSet const&);
         virtual ~TNTClusterer() {};

         void init();
         virtual void findClusters(BkgClusterCollection& clusters, ComboHitCollection const& shcol);

     private:

         unsigned formClusters(const ComboHitCollection& chcol, std::vector<BkgClusterHit>& chits, std::list<ClusterStraw>& clusters);
         void     algo1(const ComboHitCollection& chcol,std::list<ClusterStraw>& clusters, std::vector<BkgClusterHit>& chits);
         void     algo2(const ComboHitCollection& chcol,std::list<ClusterStraw>& clusters, std::vector<BkgClusterHit>& chits);

         void     initClu(const ComboHitCollection& shcol, std::vector<BkgClusterHit>& chits); 
         void     initCluMerge(const ComboHitCollection& shcol, std::vector<BkgClusterHit>& chits, std::list<ClusterStraw>& clusters); 
         void     mergeClusters(std::list<ClusterStraw>& clusters, const ComboHitCollection& chcol, float dt, float dd2);
         void     mergeTwoClu(ClusterStraw& clu1, ClusterStraw& clu2);
         float    distance(const ClusterStraw& cluster, const ComboHit& hit) const;
         
         void     dump(std::list<ClusterStraw> clusters);
         void     fillHitTree(const ComboHitCollection& chcol);
         void     fillCluTree(const ComboHitCollection& chcol, std::list<ClusterStraw>& clusters, int npass, float odist, float tdist, int nChanged);

         int          _diag;
	 bool	      _testflag;    // test background flag
         StrawHitFlag _bkgmask;     // mask for background hits
         StrawHitFlag _sigmask;     // mask for selecting signals
         bool         _comboInit;   // Start with stereo hits
         bool         _mergeInit;   // Start with stereo merge 
         float        _dseed;       // Minimum separation to seed a new cluster
         float        _dhit;        // Maximum separation to include a hit in a cluster
         float        _dd;          // cluster diameter
         float        _dt;          // natural time spread
         float        _maxdt;       // maximum time difference
         float        _maxdsum; 
         unsigned     _maxNiter;    
         unsigned     _maxNchanged;   

         std::vector<std::vector<ClusterStraw*>> _hitIndex;
         std::vector<ClusterStraw*>              _cptrs;
         float      _trms2inv; 
	 int        _ditime;
         float      _dd2; 
         float      _maxwt; 
         float      _md2;

        
         TTree* _idiag;          
         int    nhits_,hitNcombo_[8192];
         float  hitRad_[8192],hitPhi_[8192],hitTime_[8192];
         int    ncluIter_,cluId_[8192],cluNpass_[8192],cluNhit_[8192];
         float  cluRad_[8192],cluPhi_[8192],cluTime_[8192];
         float  cluR2diff_[8192],cluRdiff_[8192],cluPdiff_[8192],cluTdiff_[8192],cluTRdiff_[8192];
         int    nhitClu_, hcIdxClu_[8192], hcIdxHit_[8192], hcNpass_[8192];
         int    niter_, nclu_[8192],nChanged_[8192];
         float  odist_[8192],tdist_[8192];

	 std::map<int,int> hmap_;
	 
         
  };
}
#endif
