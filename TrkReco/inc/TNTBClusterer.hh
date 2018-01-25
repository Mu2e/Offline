//
// Two Niveau algorithm, optimized version of TLT 
//
//  Bertrand Echenard (2017) CIT
//
#ifndef TNTBClusterer_HH
#define TNTBClusterer_HH

#include "TrkReco/inc/BkgClusterer.hh"
#include "RecoDataProducts/inc/XYZVec.hh"
#include "fhiclcpp/ParameterSet.h"
#include "TTree.h"


namespace mu2e {

   struct ClusterStrawHit2 
   {      
       ClusterStrawHit2(size_t index, const ComboHit& ch, double psig2inv) : 
         _index(index),_dist(10000),_time(ch.time()),_itime(int(_time/10.0)),_pos(ch.pos()),_wdir(ch.wdir()),
         _phi(_pos.phi()),_psig2inv(psig2inv),_posResInv(1.0/ch.posRes(ComboHit::wire)),_nchanged(0)
       {}   

       size_t            _index;
       double            _dist;
       double            _time;
       int               _itime;
       XYZVec _pos;
       XYZVec _wdir;
       double            _phi;
       double            _psig2inv;
       double            _posResInv;
       unsigned          _nchanged;
   };


  class ClusterStraw2 {  
  
     public:

        ClusterStraw2(ClusterStrawHit2& hit);
        virtual ~ClusterStraw2() {};

        inline const  XYZVec& pos()              const { return _pos; }
        inline double time()                                const { return _time; }
        inline int    itime()                               const { return _itime; }
        inline bool   hasChanged()                          const { return _hasChanged;}
        inline void   flagChanged(bool val)                       { _hasChanged = val;}
        inline const  std::vector<ClusterStrawHit2*>& hits() const { return _hitsPtr; }
        inline        std::vector<ClusterStrawHit2*>& hits()       { return _hitsPtr; }
                
        void updateCache(double _maxwt);
  
  
     private:
        XYZVec             _pos;
        double                        _time;
        int                           _itime;
        std::vector<ClusterStrawHit2*> _hitsPtr; 
        bool                          _hasChanged;
  };


    
  class TNTBClusterer : public BkgClusterer
  {
     
     public:

         explicit TNTBClusterer(fhicl::ParameterSet const&);
         virtual ~TNTBClusterer() {};

         void init();
         virtual void findClusters(BkgClusterCollection& clusterColl,const ComboHitCollection& shcol);

     private:

         unsigned formClusters(std::vector<ClusterStrawHit2>& chits, std::list<ClusterStraw2>& clusters, bool init);
         double   distance(const ClusterStraw2&, ClusterStrawHit2&) const;
         double   distance2(const ClusterStraw2&, ClusterStrawHit2&) const;
         void     dump(std::list<ClusterStraw2>& clusters);

         int          _diag;
         StrawHitFlag _bkgmask; // mask for background hits
         StrawHitFlag _sigmask; // mask for selecting signals
         bool         _stereoInit;  // Start with stereo hits
         double       _dseed;       // Minimum separation to seed a new cluster
         double       _dhit;        // Maximum separation to include a hit in a cluster
         double       _dmerge;      // distance to merge 2 clusters
         double       _dd;          // cluster diameter
         double       _dt;          // natural time spread
         double       _maxdt;       // maximum time difference
         double       _trms;        // time RMS
         double       _maxdsum; 
         double       _minerr; // Spatial RMS for stereo, non-stereo hits
         double       _maxdist; // Maximum transverse distance (mm)
         double       _trms2inv, _srms2inv, _nsrms2inv; // squares of rms
         unsigned     _maxniter;    // maximum number of iterations
         unsigned     _maxnchanged;   

         std::vector<std::vector<ClusterStraw2*>> _hitIndex;
         std::vector<ClusterStraw2*> _cptrs;
	 int        _ditime;
         unsigned   _niter;  
         unsigned   _nchanged; 
         double     _md2;
         double     _dd2; 
         double     _maxwt; 
         double     _tdist;
         double     _odist;
         int        _nclu;
         int        _nhits, _nchits;
        
         TTree*     _idiag; 
         
         
  };
}
#endif
