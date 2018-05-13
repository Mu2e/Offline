//
// StrawClusterSequencePair
//
// $Id: StrawClusterSequencePair.cc,v 1.1 2013/12/07 19:51:42 brownd Exp $
// $Author: brownd $
// $Date: 2013/12/07 19:51:42 $
//
// Original author David Brown, LBNL
//
// mu2e includes
#include "TrackerMC/inc/StrawClusterSequencePair.hh"
#include "cetlib_except/exception.h"

namespace mu2e {
  namespace TrackerMC {
    StrawClusterSequencePair::StrawClusterSequencePair() {}
    StrawClusterSequencePair::StrawClusterSequencePair(StrawId sid) : 
      _scseq{StrawClusterSequence(sid,StrawEnd::cal),StrawClusterSequence(sid,StrawEnd::hv)}
    {}

    StrawClusterSequencePair::StrawClusterSequencePair(StrawClusterSequencePair const& other)
    {
      _scseq[StrawEnd::cal] = other._scseq[StrawEnd::cal];
      _scseq[StrawEnd::hv] = other._scseq[StrawEnd::hv];
    }

    StrawClusterSequencePair& StrawClusterSequencePair::operator =(StrawClusterSequencePair const& other) {
      if(&other != this){
	_scseq[StrawEnd::cal] = other._scseq[StrawEnd::cal];
	_scseq[StrawEnd::hv] = other._scseq[StrawEnd::hv];
      }
      return *this;
    }

    void StrawClusterSequencePair::insert(StrawClusterPair const& hpair) {
      if(hpair[StrawEnd::cal].strawEnd() != StrawEnd::cal || 
	  hpair[StrawEnd::hv].strawEnd() != StrawEnd::hv ||
	  hpair[0].strawId() != hpair[1].strawId())	

	throw cet::exception("SIM") 
	  << "mu2e::StrawClusterSequence: tried to add inconsistent clust pair to sequence";

      // maybe could use move symanatics here?  FIXME!!
      _scseq[StrawEnd::cal].insert(hpair[StrawEnd::cal]); 
      _scseq[StrawEnd::hv].insert(hpair[StrawEnd::hv]);
    }
  }  
}

