#ifndef SORT_FUNCTORS_HH
#define SORT_FUNCTORS_HH
//
// $Id: sort_functors.hh,v 1.2 2010/11/15 17:16:53 genser Exp $
// $Author: genser $ 
// $Date: 2010/11/15 17:16:53 $
//
// Original author KLG
//

#include "CalorimeterGeom/inc/Calorimeter.hh"

namespace mu2e {

  // utility functor to sort hits by time

  template <typename HitT>
  class lessByTime {

  public:
    
    bool operator() (HitT const & a, HitT const & b) const {
      return ( a.time() < b.time() );
    }

  };

  // utility functor to sort hits by id and time

  template <typename HitT>
  class lessByIdAndTime {

  public:
    
    bool operator() (HitT const & a, HitT const & b) const {
      return (a.id() < b.id() ||
              (a.id() == b.id() &&
               a.time() < b.time() 
               ) 
              );
    }

  };

  // utility functor to sort hits by crystal id and time

  template <typename HitT>
  class lessByCIdAndTime {

  public:
    

    explicit lessByCIdAndTime(Calorimeter const & cal): _cal(cal) {}

    bool operator() (HitT const & a, HitT const & b) const {
      
      return ( _cal.getCrystalByRO(a.id()) < _cal.getCrystalByRO(b.id()) ||
               (_cal.getCrystalByRO(a.id()) == _cal.getCrystalByRO(b.id()) &&
                a.time() < b.time() 
                ) 
               );
    }

  private:

    Calorimeter const & _cal;

  };
  
}

#endif
