#ifndef TrackerMC_StrawClusterSequencePair_hh
#define TrackerMC_StrawClusterSequencePair_hh
// Define ClusterSequencePair, ie the clust sequence from both ends of a straw
//
// $Id: StrawClusterSequencePair.hh,v 1.1 2013/12/07 19:51:42 brownd Exp $
// $Author: brownd $
// $Date: 2013/12/07 19:51:42 $
//
// Original author David Brown, LBNL
#include "TrackerMC/inc/StrawClusterSequence.hh"
#include <utility>
namespace mu2e {
  namespace TrackerMC {
    class StrawClusterSequencePair{
      public:
	typedef StrawCluster StrawClusterPair[2];
	StrawClusterSequencePair();
	StrawClusterSequencePair(StrawIndex index);
	StrawClusterSequencePair(StrawClusterSequencePair const& other);
	StrawClusterSequencePair& operator =(StrawClusterSequencePair const& other);
	StrawClusterSequence& clustSequence(StrawEnd end) { if(end==StrawEnd::plus) return _plus; else return _minus; }
	StrawClusterSequence const& clustSequence(StrawEnd end) const { if(end==StrawEnd::plus) return _plus; else return _minus; }
	void insert(StrawClusterPair const& hpair);
	StrawIndex strawIndex() const { return _plus.strawIndex(); }
      private:
	StrawClusterSequence _plus, _minus; // plus and minus end sequences
    };
  }
}
#endif
