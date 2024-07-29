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
        return uRes();
      }
      case ComboHit::trans : {
        return vRes();
      }
      case ComboHit::z : {
        return wRes();
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

  ComboHitCollection::CHCPTR ComboHitCollection::parent(StrawIdMask::Level level, bool stopatfirst) const {
    auto retval = _parent;
    if(_parent.refCore().isNull()){
      throw cet::exception("RECO")<<"mu2e::ComboHitCollection: no such parent" << std::endl;
    } else {
      // recursive call
      if (retval->level() != level || (!stopatfirst && retval->parent().refCore().isNonnull() && retval->parent()->level() == level)){
        retval = _parent->parent(level);
      }
    }
    return retval;
  }

  void ComboHitCollection::setAsSubset(art::Handle<ComboHitCollection> const& ohandle){
    if (ohandle.isValid()){
      auto optr = CHCPTR(ohandle);
      setAsSubset(optr);
    } else {
      throw cet::exception("RECO")<<"mu2e::ComboHitCollection: invalid handle" << std::endl;
    }
  }
  void ComboHitCollection::setAsSubset(art::ValidHandle<ComboHitCollection> const& ohandle){
    auto optr = CHCPTR(ohandle);
    setAsSubset(optr);
  }
  void ComboHitCollection::setAsSubset(CHCPTR const& other){
    _parent = other->parent();
    if (_parent.refCore().isNull()){
      _parent = other;
    }
  }

  void ComboHitCollection::setAsSubset(art::ValidHandle<ComboHitCollection> const& ohandle, StrawIdMask::Level level, bool stopatfirst){
    if (ohandle.isValid()){
      auto optr = CHCPTR(ohandle);
      setAsSubset(optr,level,stopatfirst);
    } else {
      throw cet::exception("RECO")<<"mu2e::ComboHitCollection: invalid handle" << std::endl;
    }
  }
  void ComboHitCollection::setAsSubset(art::Handle<ComboHitCollection> const& ohandle, StrawIdMask::Level level, bool stopatfirst){
    auto optr = CHCPTR(ohandle);
    setAsSubset(optr,level,stopatfirst);
  }
  void ComboHitCollection::setAsSubset(CHCPTR const& optr, StrawIdMask::Level level, bool stopatfirst){
    if (optr->level() == level && (stopatfirst || optr->parent().refCore().isNull() || optr->parent()->level() != level)){
      _parent = optr;
    }else{
      _parent = optr->parent(level,stopatfirst);
    }
  }


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

  void ComboHitCollection::fillStrawDigiIndices(size_t chindex, vector<StrawDigiIndex>& shiv, bool stopatfirst) const {
    // if this is a straw-level ComboHit, get the digi index
    ComboHit const& ch = this->at(chindex);
    if(level() == StrawIdMask::uniquestraw && (stopatfirst || _parent.refCore().isNull() || _parent->level() != StrawIdMask::uniquestraw)){
      shiv.push_back(ch.index(0));
      // if this collection references a sub-collection, go down recursively
    } else if(_parent.refCore().isNonnull()){
      for(size_t iind = 0;iind < ch.nCombo(); ++iind){
        _parent->fillStrawDigiIndices(ch.index(iind),shiv);
      }
    } else {
      throw cet::exception("RECO")<<"mu2e::ComboHitCollection: invalid ComboHit" << std::endl;
    }
  }

  ComboHitCollection const* ComboHitCollection::fillStrawHitIndices( size_t chindex, SHIV& shiv, StrawIdMask::Level clevel, bool stopatfirst) const {
    ComboHitCollection const* retval = this;
    if(level() == clevel && (stopatfirst || _parent.refCore().isNull() || _parent->level() != clevel)){
      shiv.push_back(chindex);
      // if this collection references other collections, go down recursively
    } else if(_parent.refCore().isNonnull()){
      ComboHit const& ch = this->at(chindex);
      if(_parent->level() == clevel){
        retval = _parent.get();
        shiv.reserve(shiv.size() + ch.nCombo());
        for(size_t iind = 0;iind < ch.nCombo(); ++iind){
          shiv.push_back(ch.index(iind));
        }
      } else {
        SHIV tempshiv;
        tempshiv.reserve(2*ch.nCombo());
        for(size_t iind = 0;iind < ch.nCombo(); ++iind){
          tempshiv.push_back(ch.index(iind));
        }
        retval = _parent->fillStrawHitIndices(tempshiv,shiv,clevel);
      }
    } else {
      throw cet::exception("RECO")<<"mu2e::ComboHitCollection: invalid ComboHitCollection nesting" << std::endl;
    }
    return retval;
  }

  ComboHitCollection const* ComboHitCollection::fillStrawHitIndices(SHIV const& inshiv, SHIV& outshiv, StrawIdMask::Level clevel, bool stopatfirst) const {
    ComboHitCollection const* retval = this;
    if(level() == clevel && (stopatfirst || _parent.refCore().isNull() || _parent->level() != clevel)){
      outshiv = inshiv;
    } else if(_parent.refCore().isNonnull()){
      if(_parent->level() == clevel){
        // parent hits are at the correct level
        retval = _parent.get();
        outshiv.reserve(2*inshiv.size()); // just a guess TODO
        for(auto shi : inshiv) {
          auto const& ch = (*this)[shi]; // don't need to check range
          for(size_t iind = 0;iind < ch.nCombo(); ++iind){
            outshiv.push_back(ch.index(iind));
          }
        }
      } else {
        SHIV tempshiv;
        tempshiv.reserve(2*inshiv.size());
        for(auto shi : inshiv) {
          auto const& ch = (*this)[shi]; // don't need to check range
          for(size_t iind = 0;iind < ch.nCombo(); ++iind){
            tempshiv.push_back(ch.index(iind));
          }
        }
        retval = _parent->fillStrawHitIndices(tempshiv,outshiv,clevel);
      }
    } else {
      throw cet::exception("RECO")<<"mu2e::ComboHitCollection: invalid ComboHitCollection nesting" << std::endl;
    }
    return retval;
  }

  ComboHitCollection const* ComboHitCollection::fillStrawHitIndices( SHIV& shiv, StrawIdMask::Level clevel, bool stopatfirst) const {
    shiv.clear();
    const ComboHitCollection* retval = this;
    if(level() == clevel && (stopatfirst || _parent.refCore().isNull() || _parent->level() != clevel)){
      shiv.reserve(this->size());
      for(size_t iind = 0;iind < this->size(); ++iind)shiv.push_back(iind);
      // if this collection references other collections, go down recursively
    } else if(_parent.refCore().isNonnull()){
      if(_parent->level() == clevel){
        retval = _parent.get();
        shiv.reserve(2*(this->size())); // estimated combo factor
        // current hits reference into the desired level.  Fill the indices from the current hits and done
        for(auto const& ch : *this){
          for(size_t iind = 0;iind < ch.nCombo(); ++iind){
            shiv.push_back(ch.index(iind));
          }
        }
      } else {
        // we need to descend, but only picking up the hits in the sub-collection referenced by my hits
        // build a temporary array
        SHIV tempshiv;
        tempshiv.reserve(2*this->size()); // size is a guess, verify on data TODO
        for(auto const& ch : *this){
          for(size_t iind = 0;iind < ch.nCombo(); ++iind){
            tempshiv.push_back(ch.index(iind));
          }
        }
        retval = _parent->fillStrawHitIndices(tempshiv,shiv,clevel);
      }
    } else {
      throw cet::exception("RECO")<<"mu2e::ComboHitCollection: invalid ComboHitCollection Nesting" << std::endl;
    }
    return retval;
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

  void ComboHitCollection::fillComboHits( size_t chindex, CHCIter& iters) const {
    iters.clear();
      // if this collection references other collections: if so, we can fill the vector
    if(_parent.refCore().isNonnull()){
      ComboHit const& ch = (*this)[chindex];
      for(size_t iind = 0;iind < ch.nCombo(); ++iind){
        iters.push_back(std::next(_parent->begin(), ch.index(iind)));
      }
    }
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
