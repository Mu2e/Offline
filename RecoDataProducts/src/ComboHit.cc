//
// Class to describe a combination of Straw Hits
//
// Original author David Brown
//
// Mu2e includes
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
// art includes
#include "cetlib_except/exception.h"
// c++ includes
#include <iostream>
#include <limits>
using std::vector;
namespace mu2e {

  Float_t ComboHit::posRes(edir dir) const {
    switch ( dir ) {
      case ComboHit::wire : {
        return _ures;
      }
      case ComboHit::trans : {
        return _vres;
      }
      case ComboHit::z : {
        return _wres;
      }
      default : {
        return -1.0;
      }
    }
  }

  void ComboHit::init(ComboHit const& other, size_t index) {
    if(index > std::numeric_limits<size_t>::max())
      throw cet::exception("RECO")<<"mu2e::ComboHitCollection: invalid index" << std::endl;
    *this = other;
    _ncombo = 1;
    _pind[0] = index;
  }

  bool ComboHit::addIndex(size_t shi) {
    if(shi < std::numeric_limits<size_t>::max() && _ncombo < MaxNCombo){
      _pind[_ncombo] = shi;
      ++_ncombo;
      return true;
    } else
      return false;
  }

#ifndef __ROOTCLING__
  void ComboHitCollection::setSameParent(ComboHitCollection const& other) { _parent = other.parent(); }
  void ComboHitCollection::setParent(CHCPTR const& parent) { _parent = parent; }

  void ComboHitCollection::setParent(art::Handle<ComboHitCollection> const& phandle) {
    if(phandle.isValid()){
      _parent = CHCPTR(phandle);
    } else {
      throw cet::exception("RECO")<<"mu2e::ComboHitCollection: invalid handle" << std::endl;
    }
  }

  void ComboHitCollection::setParent(art::ValidHandle<ComboHitCollection> const& phandle) {
    _parent = CHCPTR(phandle);
  }

  void ComboHitCollection::fillStrawDigiIndices(size_t chindex, vector<StrawDigiIndex>& shids) const {
    // if this is a straw-level ComboHit, get the digi index
    ComboHit const& ch = this->at(chindex);
    if(level() == StrawIdMask::uniquestraw) {
      shids.push_back(ch.index(0));
      // if this collection references a sub-collection, go down recursively
    } else if(_parent.refCore().isNonnull()){
      for(size_t iind = 0;iind < ch.nCombo(); ++iind){
        _parent->fillStrawDigiIndices(ch.index(iind),shids);
      }
    } else {
      throw cet::exception("RECO")<<"mu2e::ComboHitCollection: invalid ComboHit" << std::endl;
    }
  }

  void ComboHitCollection::fillStrawHitIndices( size_t chindex, vector<StrawHitIndex>& shids) const {
    // if this is a straw-level Collection, the input is the index
    ComboHit const& ch = this->at(chindex);
    if(level() == StrawIdMask::uniquestraw) {
      shids.push_back(chindex);
      // if this collection references other collections, go down recursively
    } else if(_parent.refCore().isNonnull()){
       for(size_t iind = 0;iind < ch.nCombo(); ++iind){
        _parent->fillStrawHitIndices(ch.index(iind),shids);
      }
    } else {
      throw cet::exception("RECO")<<"mu2e::ComboHitCollection: invalid ComboHit" << std::endl;
    }
  }

  void ComboHitCollection::fillStrawHitIndices( vector<vector<StrawHitIndex> >& shids) const {
    // reset
    shids = vector<vector<StrawHitIndex> >(size());
    // if this is a straw-level Collection, the input is the index
    if(level() == StrawIdMask::uniquestraw) {
      for(size_t ich = 0;ich < size(); ++ich){
        shids[ich].push_back(ich);
      }
    } else if(_parent.refCore().isNonnull()){
      vector<vector<StrawHitIndex>> pshids(size());
      _parent->fillStrawHitIndices(pshids);
        // check down 1 more layer
      if(_parent->parent().refCore().isNonnull()){
        // call down 1 layer
        vector<vector<StrawHitIndex>> pshids(size());
        _parent->fillStrawHitIndices(pshids);
        // roll these up
        for(size_t ich=0;ich < size();++ich){
          ComboHit const& ch = (*this)[ich];
          for(auto iph : ch.indexArray())
            shids[ich].insert(shids[ich].end(),pshids[iph].begin(),pshids[iph].end());//FIX ME , problem with indices when using StereoHits, side stepped by if(pshids.size() == 8)
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
  }

  void ComboHitCollection::fillComboHits( std::vector<uint16_t> const& indices, CHCIter& iters) const {
    if(_parent.refCore().isNull()){
      for(auto index : indices){
        iters.push_back(std::next(begin(),index));
      }
    } else {
      std::vector<uint16_t> subindices;
      for(auto index : indices){
        ComboHit const& ch = (*this)[index];
        for(size_t iind = 0;iind < ch.nCombo(); ++iind)
          subindices.push_back(ch.index(iind));
        }
        _parent->fillComboHits(subindices,iters);
    }
  }

  bool ComboHitCollection::fillComboHits( size_t chindex, CHCIter& iters) const {
    bool retval(false);
    iters.clear();
      // if this collection references other collections: if so, we can fill the vector
    if(_parent.refCore().isNonnull()){
      ComboHit const& ch = (*this)[chindex];
      for(size_t iind = 0;iind < ch.nCombo(); ++iind){
        iters.push_back(std::next(_parent->begin(), ch.index(iind)));
      }
      retval = true;
    }
    return retval;
  }

#endif
  StrawIdMask::Level ComboHitCollection::level() const { return this->size() > 0 ? this->front().mask().level() : StrawIdMask::none; }

  unsigned ComboHitCollection::nStrawHits() const {
    unsigned retval(0);
    for(auto const& ch : *this )
      retval += ch.nStrawHits();
    return retval;
  }

  void ComboHit::print( std::ostream& ost, bool doEndl) const {
    ost << " ComboHit:"
        << " id "      << _sid
        << " level "   << _mask.level()
        << " time "     << _time
        << " wdist " << _wdist
        << " position " << _pos
        << " early end " << _eend
        << " flag " << _flag
        << " edep "     << _edep
        << " ncombo " << _ncombo
        << " nStrawHit " << _nsh;

    if ( doEndl ){
      ost << std::endl;
    }

  }
}
