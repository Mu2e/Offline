#ifndef SORT_FUNCTORS_HH
#define SORT_FUNCTORS_HH
//
// $Id: sort_functors.hh,v 1.1 2010/11/12 19:40:51 genser Exp $
// $Author: genser $ 
// $Date: 2010/11/12 19:40:51 $
//
// Original author KLG
//

namespace mu2e {

  // utility functor to sort hits by time

  template <typename HitT>
  class lessByTime {

  public:
    
    bool operator() (const HitT& a, const HitT& b) const {
      return ( a.time() < b.time() );
    }

  };

  // utility functor to sort hits by id & time

  template <typename HitT>
  class lessByIdAndTime {

  public:
    
    bool operator() (const HitT& a, const HitT& b) const {
      return (a.id() < b.id() ||
              (a.id() == b.id() &&
               a.time() < b.time() 
               ) 
              );
    }

  };

  // utility functor to sort hits by crystal id & time

  template <typename HitT>
  class lessByCIdAndTime {

  public:
    
    // explicit lessByCIdAndTime(GeomHandle<Calorimeter> cg) _cg(cg) {} 
    // GeomHandle is not copyable

    bool operator() (const HitT& a, const HitT& b) const {
      
      GeomHandle<Calorimeter> _cg;
      return ( _cg->getCrystalByRO(a.id()) < _cg->getCrystalByRO(b.id()) ||
               (_cg->getCrystalByRO(a.id()) == _cg->getCrystalByRO(b.id()) &&
                a.time() < b.time() 
                ) 
               );
    }

  };
  
}

#endif
