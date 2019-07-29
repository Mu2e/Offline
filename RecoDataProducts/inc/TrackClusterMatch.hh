#ifndef TrackCaloMatching_TrackClusterMatch_hh
#define TrackCaloMatching_TrackClusterMatch_hh

// Original author Gianantonio Pezzullo
#include <ostream>

#include "canvas/Persistency/Common/Ptr.h"
#include "RecoDataProducts/inc/TrkCaloIntersect.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"

namespace mu2e {

  typedef art::Ptr<TrkCaloIntersect> TrkCaloIntersectPtr ;
  typedef art::Ptr<CaloCluster     > CaloClusterPtr ;

  class TrackClusterMatch {
  public:

    struct Data_t {
      int      icl;			// cluster index
      int      iex;			// extrapolated track index
      double   xtrk;			// track coordinates
      double   ytrk;
      double   ztrk;
      double   ttrk;			// track time
      double   nx;			// track direction in the "interaction point"
      double   ny;
      double   nz;
      double   dx;
      double   dy;
      double   dz;
      double   du;
      double   dv;
      double   dt;
      double   ep;
      double   chi2;
      double   chi2_time;
      double   int_depth;
      double   ds;			// path length inside the disk
      double   dr;                      // R(cluster)-R(track)
      double   sint;                    // "interaction length"
    };

  protected:
    int                   _icl;		// cluster index
    int                   _iex;		// extrapolated track index
    TrkCaloIntersectPtr   _textrapol; 
    CaloClusterPtr        _cluster; 
    double                _xtrk;	// track coordinates
    double                _ytrk;
    double                _ztrk;
    double                _ttrk;
    double                _nx;		// track direction
    double                _ny;
    double                _nz;
    double                _dx;
    double                _dy;
    double                _dz;
    double                _du;
    double                _dv;
    double                _dt;
    double                _ep;
    double                _chi2;
    double                _chi2_time;
    double                _int_depth;
    double                _ds;
    double                _dr;
    double                _sint;

  public:
    TrackClusterMatch();
    TrackClusterMatch(TrkCaloIntersectPtr& Textrapol, CaloClusterPtr & Cluster, Data_t* Data);
    ~TrackClusterMatch();

    int                       icl        () const { return _icl; }
    int                       iex        () const { return _iex; }
    const TrkCaloIntersect*   textrapol  () const { return _textrapol.get(); }
    const CaloCluster*        caloCluster() const { return _cluster.get(); }
    double                    xtrk       () const { return _xtrk; }
    double                    ytrk       () const { return _ytrk; }
    double                    ztrk       () const { return _ztrk; }
    double                    ttrk       () const { return _ttrk; }
    double                    nx         () const { return _nx; }
    double                    ny         () const { return _ny; }
    double                    nz         () const { return _nz; }
    double                    dx         () const { return _dx; }
    double                    dy         () const { return _dy; }
    double                    dz         () const { return _dz; }
    double                    dt         () const { return _dt; }
    double                    du         () const { return _du; }
    double                    dv         () const { return _dv; }
    double                    ep         () const { return _ep; }
    double                    chi2       () const { return _chi2; }
    double                    chi2_time  () const { return _chi2_time; }
    double                    int_depth  () const { return _int_depth; }
    double                    ds         () const { return _ds; }
    double                    dr         () const { return _dr; }
    double                    sint       () const { return _sint; }

    void print(const char* Option) const ;
  };

  typedef std::vector<mu2e::TrackClusterMatch> TrackClusterMatchCollection;

}


#endif/*TrackCaloMatching_TrackClusterMatch_hh*/


