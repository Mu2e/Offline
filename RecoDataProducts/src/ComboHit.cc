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

  ComboHit::ComboHit() : _wres(-1.0),_tres(-1.0), _wdist(0.), _time(0.0), _edep(0.0), _qual(0.0), _ncombo(0), _nsh(0) {}

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
    std::vector<art::Handle<ComboHitCollection> > all_handles;
    // exhaustive search is fast enough
    event.getManyByType(all_handles);
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

  void ComboHitCollection::fillStrawIds(art::Event const& event, uint16_t chindex, std::vector<StrawHitIndex>& shids) const {
    ComboHit const& ch = (*this)[chindex];
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
	  pc->fillStrawIds(event,ch.index(iind),shids);
	}
      } else {
	throw cet::exception("RECO")<<"mu2e::ComboHitCollection: Can't find parent collection" << std::endl;
      }
    } else {
      // if not, it is the bottom and references StrawHits; fill the index vector with the content
      for(uint16_t iind = 0;iind < ch.nCombo(); ++iind){
	shids.push_back(ch.index(iind));
      }
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
}
