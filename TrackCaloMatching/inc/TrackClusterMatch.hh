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
  protected:
    TrkToCaloExtrapolPtr  _textrapol; 
    CaloClusterPtr        _cluster; 
    double                _chi2;

  public:
    TrackClusterMatch();
    TrackClusterMatch(TrkToCaloExtrapolPtr& Textrapol, CaloClusterPtr & Cluster, double Chi2);
    ~TrackClusterMatch();

    const TrkToCaloExtrapol*  textrapol  () const { return _textrapol.get(); }
    const CaloCluster*        caloCluster() const { return _cluster.get(); }
    double                    chi2       () const { return _chi2; }

    void print(const char* Option) const ;
  };

  typedef std::vector<mu2e::TrackClusterMatch> TrackClusterMatchCollection;

}


#endif/*TrackCaloMatching_TrackClusterMatch_hh*/


