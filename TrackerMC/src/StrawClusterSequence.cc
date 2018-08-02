//
// StrawClusterSequence is a time-ordered sequence of StrawClusters
//
// Original author David Brown, LBNL
//
// mu2e includes
#include "TrackerMC/inc/StrawClusterSequence.hh"
#include "cetlib_except/exception.h"

using namespace std;

namespace mu2e {
  namespace TrackerMC {
    StrawClusterSequence::StrawClusterSequence() : _strawId(0), _end(StrawEnd::cal)
    {}

    StrawClusterSequence::StrawClusterSequence(StrawCluster const& clust) :
      _strawId(clust.strawId()), _end(clust.strawEnd())
    {
      insert(clust);
    }

    StrawClusterSequence::StrawClusterSequence(StrawId const& sid, StrawEnd end) :
      _strawId(sid), _end(end)
    {}

    StrawClusterSequence::StrawClusterSequence(StrawClusterSequence const& other) :
      _strawId(other._strawId), _end(other._end), _clist(other._clist) {}

    StrawClusterSequence& StrawClusterSequence::operator =(StrawClusterSequence const& other) {
      if(&other != this){
	_strawId = other._strawId;
	_end = other._end;
	_clist = other._clist;
      }
      return *this;
    }
    // insert a new clust.  This is the only non-trivial function
    ClusterList::iterator StrawClusterSequence::insert(StrawCluster const& clust) {
      ClusterList::iterator retval = _clist.end();
      if(clust.type() == StrawCluster::unknown){
	throw cet::exception("SIM")
	  << "mu2e::StrawClusterSequence: tried to add unknown clust type"
	  << endl;
	return retval;
      }
      // make sure the straw and end are the same
      if(!_clist.empty() && (clust.strawId() != strawId()
	    || clust.strawEnd() != strawEnd())){
	throw cet::exception("SIM")
	  << "mu2e::StrawClusterSequence: tried to add clust from a different straw/end to a sequence"
	  << endl;
	return retval;
      }
      if(_clist.empty()){
	_clist.push_front(clust);
	_strawId = clust.strawId();
	_end = clust.strawEnd();
	retval = _clist.begin();
      } else {
	// loop over the contents and insert in the correct place
	ClusterList::iterator ibefore = _clist.begin();
	while(ibefore != _clist.end() && ibefore->time() < clust.time())
	  ++ibefore;
	retval = _clist.insert(ibefore,clust);
      }
      return retval;
    }
  }
}
