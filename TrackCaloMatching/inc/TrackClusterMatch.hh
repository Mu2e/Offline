#ifndef TrackCaloMatching_TrackClusterMatch_hh
#define TrackCaloMatching_TrackClusterMatch_hh

// Original author Gianantonio Pezzullo
#include <ostream>

#include "art/Persistency/Common/Ptr.h"
#include "TrackCaloMatching/inc/TrkToCaloExtrapol.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"

namespace mu2e {

  typedef art::Ptr<TrkToCaloExtrapol> TrkToCaloExtrapolPtr ;
  typedef art::Ptr<CaloCluster      > CaloClusterPtr ;

  class TrackClusterMatch {
  public:

    struct Data_t {
      int      icl;			// cluster index
      int      iex;			// extrapolated track index
      double   dx;
      double   dy;
      double   dz;
      double   du;
      double   dv;
      double   dt;
      double   ep;
      double   chi2;
    };

  protected:
    TrkToCaloExtrapolPtr  _textrapol; 
    CaloClusterPtr        _cluster; 
    double                _dx;
    double                _dy;
    double                _dz;
    double                _du;
    double                _dv;
    double                _dt;
    double                _ep;
    double                _chi2;

  public:
    TrackClusterMatch();
    TrackClusterMatch(TrkToCaloExtrapolPtr& Textrapol, CaloClusterPtr & Cluster, Data_t* Data);
    ~TrackClusterMatch();

    const TrkToCaloExtrapol*  textrapol  () const { return _textrapol.get(); }
    const CaloCluster*        caloCluster() const { return _cluster.get(); }
    double                    dx         () const { return _dx; }
    double                    dy         () const { return _dy; }
    double                    dz         () const { return _dz; }
    double                    dt         () const { return _dt; }
    double                    du         () const { return _du; }
    double                    dv         () const { return _dv; }
    double                    ep         () const { return _ep; }
    double                    chi2       () const { return _chi2; }

    void print(const char* Option) const ;
  };

  typedef std::vector<mu2e::TrackClusterMatch> TrackClusterMatchCollection;

}


#endif/*TrackCaloMatching_TrackClusterMatch_hh*/


