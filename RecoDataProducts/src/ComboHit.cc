//
// Class to describe a combination of Straw Hits
//
// Original author David Brown
//
// Mu2e includes
#include "RecoDataProducts/inc/ComboHit.hh"
// art includes
#include "cetlib_except/exception.h"
// c++ includes
#include <iostream>
using std::vector;
namespace mu2e {

  Float_t ComboHit::posRes(edir dir) const {
    switch ( dir ) {
      case ComboHit::wire : {
	return _wres;
      }
      case ComboHit::trans : {
	return _tres;
      }
      default : {
	return -1.0;
      }
    }
  }

  ComboHit::ComboHit() : _wres(-1.0),_tres(-1.0), _wdist(0.), _time(0.0), _edep(0.0), _qual(0.0), _dtime(0.0), _ptime(0.0), _pathlength(0.0), _hphi(0.0), _ncombo(0), _nsh(0), _pind{0} {}

  ComboHit::ComboHit(ComboHit const& shp, StrawDigiIndex hi, double phi) : ComboHit(shp)
  {
    _hphi = phi;
  }

  void ComboHit::init(ComboHit const& other, uint16_t index) {
    *this = other;
    _ncombo = 1;
    _pind[0] = index;
  }

  uint16_t ComboHit::index(uint16_t ish) const {
    if(ish < _ncombo)
      return _pind[ish];
    else
      throw cet::exception("RECO")<<"mu2e::ComboHit: invalid index" << std::endl;
  }

  bool ComboHit::addIndex(uint16_t shi) {
    if(_ncombo < MaxNCombo){
      _pind[_ncombo] = shi;
      ++_ncombo;
      return true;
    } else
      return false;
  }

  void ComboHitCollection::setParentHandle(art::Event const& event, art::Handle<ComboHitCollection>& phandle) const  {
    // set the handle to an invalid state in case we find no such
    phandle = art::Handle<ComboHitCollection>();
    vector<art::Handle<ComboHitCollection> > all_handles =  event.getMany<ComboHitCollection>();
    // exhaustive search is fast enough
    for (auto const& handle : all_handles) {
      if(_parent == handle.id()){
	phandle = handle;
	break;
      }
    }
  }

  void ComboHitCollection::setParent(art::Handle<ComboHitCollection> const& phandle) {
    if(phandle.isValid()){
      _parent = phandle.id();
    } else {
      throw cet::exception("RECO")<<"mu2e::ComboHitCollection: invalid handle" << std::endl;
    }
  }

  void ComboHitCollection::fillStrawDigiIndices(art::Event const& event, uint16_t chindex, vector<StrawDigiIndex>& shids) const {
    ComboHit const& ch = this->at(chindex);
   // see if this collection references other collections: if so, go down 1 layer
    if(_parent.isValid()){
    // get the parent handle
      art::Handle<ComboHitCollection> ph;
      setParentHandle(event,ph);
      if(ph.isValid()){
      // get the parent collection
	const ComboHitCollection *pc = ph.product();
	// recursive calls on the ComboHits in the parent collection referenced by the specified ComboHit
	for(uint16_t iind = 0;iind < ch.nCombo(); ++iind){
	  pc->fillStrawDigiIndices(event,ch.index(iind),shids);
	}
      } else {
	throw cet::exception("RECO")<<"mu2e::ComboHitCollection: Can't find parent collection" << std::endl;
      }
    } else {
      if(ch.nCombo() != 1 || ch.nStrawHits() != 1)
	throw cet::exception("RECO")<<"mu2e::ComboHitCollection: invalid ComboHit" << std::endl;
      // if not, it is the bottom and references StrawHits; fill the index vector with the content
      shids.push_back(ch.index(0));
    }
  }

  void ComboHitCollection::fillStrawHitIndices(art::Event const& event, uint16_t chindex, vector<StrawHitIndex>& shids) const {
    ComboHit const& ch = this->at(chindex);
   // see if this collection references other collections: if so, go down 1 layer
    if(_parent.isValid()){
    // get the parent handle
      art::Handle<ComboHitCollection> ph;
      setParentHandle(event,ph);
      if(ph.isValid()){
      // get the parent collection
	const ComboHitCollection *pc = ph.product();
	// recursive calls on the ComboHits in the parent collection referenced by the specified ComboHit
	for(uint16_t iind = 0;iind < ch.nCombo(); ++iind){
	  pc->fillStrawHitIndices(event,ch.index(iind),shids);
	}
      } else {
	throw cet::exception("RECO")<<"mu2e::ComboHitCollection: Can't find parent collection" << std::endl;
      }
    } else {
      if(ch.nCombo() != 1 || ch.nStrawHits() != 1)
	throw cet::exception("RECO")<<"mu2e::ComboHitCollection: invalid ComboHit" << std::endl;
      // if not, it is the bottom and the combo hit index is the same as the StrawHit index
      shids.push_back(chindex);
    }
  }

  void ComboHitCollection::fillStrawHitIndices(art::Event const& event, vector<vector<StrawHitIndex> >& shids) const {
  // reset
    shids = vector<vector<StrawHitIndex> >(size());
   // see if this collection references other collections: if so, go down 1 layer
    if(parent().isValid()){
    // get the parent handle
      art::Handle<ComboHitCollection> ph;
      setParentHandle(event,ph);
      if(ph.isValid()){
	const ComboHitCollection *pc = ph.product();
	// check down 1 more layer
	if(pc->parent().isValid()){
	// call down 1 layer
	  vector<vector<StrawHitIndex>> pshids(size());
	  pc->fillStrawHitIndices(event,pshids);
	// roll these up
	  for(size_t ich=0;ich < size();++ich){
	    ComboHit const& ch = (*this)[ich];
	    for(auto iph : ch.indexArray())
	      shids[ich].insert(shids[ich].end(),pshids[iph].begin(),pshids[iph].end());
	  }
	} else {
	  // can skip a step in the hierarchy since the parent of this collection is at the top
	  for(size_t ich=0;ich < size();++ich){
	    ComboHit const& ch = (*this)[ich];
	    shids[ich].insert(shids[ich].end(),ch.indexArray().begin(),ch.indexArray().end());
	  }
	}
      } else {
	throw cet::exception("RECO")<<"mu2e::ComboHitCollection: Can't find parent collection" << std::endl;
      }
    } else {
    // already at the top: fill with self index
      for(size_t ich = 0;ich < size(); ++ich)
	shids[ich].push_back(ich);
    }
  }

  void ComboHitCollection::fillComboHits(art::Event const& event, std::vector<uint16_t> const& indices, CHCIter& iters) const {
    if(_parent.isValid()){
    // get the parent handle
      art::Handle<ComboHitCollection> ph;
      setParentHandle(event,ph);
      if(ph.isValid()){
      // get the parent collection
	const ComboHitCollection *pc = ph.product();
	// translate the indices down
	std::vector<uint16_t> subindices;
	for(auto index : indices){
	  ComboHit const& ch = (*this)[index];
	  for(uint16_t iind = 0;iind < ch.nCombo(); ++iind)
	    subindices.push_back(ch.index(iind));
	}
	pc->fillComboHits(event,subindices,iters);
      } else {
	throw cet::exception("RECO")<<"mu2e::ComboHitCollection: Can't find parent collection" << std::endl;
      }
    } else {
// already at lowest level: just translate indices to pointers
      for(auto index : indices)
	iters.push_back(std::next(begin(),index));
    }
  }
 
  bool ComboHitCollection::fillComboHits(art::Event const& event, uint16_t chindex, CHCIter& iters) const {
    bool retval(false);
    iters.clear();
    ComboHit const& ch = (*this)[chindex];
   // see if this collection references other collections: if so, we can fill the vector, if not reference myself
    if(_parent.isValid()){
    // get the parent handle
      art::Handle<ComboHitCollection> ph;
      setParentHandle(event,ph);
      if(ph.isValid()){
      // get the parent collection
	const ComboHitCollection *pc = ph.product();
	for(uint16_t iind = 0;iind < ch.nCombo(); ++iind){
	  iters.push_back(std::next(pc->begin(), ch.index(iind)));
	}
	retval = true;
      } else {
	throw cet::exception("RECO")<<"mu2e::ComboHitCollection: Can't find parent collection" << std::endl;
      }
    }
    return retval; 
  }

  uint16_t ComboHitCollection::nStrawHits() const {
    uint16_t retval(0);
    for(auto const& ch : *this )
      retval += ch.nStrawHits();
    return retval;
  }

  void ComboHit::print( std::ostream& ost, bool doEndl) const {
    ost << " ComboHit:"
        << " id "      << _sid
        << " time "     << _time
        << " drift time " << _dtime
        << " prop time " << _ptime
        << " path length " << _pathlength
	<< " position " << _pos
	<< " end " << _tend
	<< " flag " << _flag
        << " eDep "     << _edep
	<< " ncombo " << _ncombo
	<< " nStrawHit " << _nsh;

    if ( doEndl ){
      ost << std::endl;
    }


  }
}

