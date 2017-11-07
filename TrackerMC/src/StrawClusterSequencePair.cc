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
    StrawClusterSequencePair::StrawClusterSequencePair(StrawIndex index) : 
      _scseq{StrawClusterSequence(index,TrkTypes::cal),StrawClusterSequence(index,TrkTypes::hv)}
    {}

    StrawClusterSequencePair::StrawClusterSequencePair(StrawClusterSequencePair const& other)
    {
      _scseq[TrkTypes::cal] = other._scseq[TrkTypes::cal];
      _scseq[TrkTypes::hv] = other._scseq[TrkTypes::hv];
    }

    StrawClusterSequencePair& StrawClusterSequencePair::operator =(StrawClusterSequencePair const& other) {
      if(&other != this){
	_scseq[TrkTypes::cal] = other._scseq[TrkTypes::cal];
	_scseq[TrkTypes::hv] = other._scseq[TrkTypes::hv];
      }
      return *this;
    }

    void StrawClusterSequencePair::insert(StrawClusterPair const& hpair) {
      if(hpair[TrkTypes::cal].strawEnd() != TrkTypes::cal || 
	  hpair[TrkTypes::hv].strawEnd() != TrkTypes::hv ||
	  hpair[0].strawIndex() != hpair[1].strawIndex())	

	throw cet::exception("SIM") 
	  << "mu2e::StrawClusterSequence: tried to add inconsistent clust pair to sequence";

      // maybe could use move symanatics here?  FIXME!!
      _scseq[TrkTypes::cal].insert(hpair[TrkTypes::cal]); 
      _scseq[TrkTypes::hv].insert(hpair[TrkTypes::hv]);
    }
  }  
}

