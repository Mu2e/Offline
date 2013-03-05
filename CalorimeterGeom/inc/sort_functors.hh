#ifndef CalorimeterGeom_sort_functors_hh
#define CalorimeterGeom_sort_functors_hh
//
// $Id: sort_functors.hh,v 1.5 2013/03/05 20:33:25 aluca Exp $
// $Author: aluca $
// $Date: 2013/03/05 20:33:25 $
//
// Original author KLG
//

#include "CalorimeterGeom/inc/Calorimeter.hh"

namespace mu2e {

  // utility functor to sort hits by time
  template <typename HitT>  class lessByTime {

    public:

      bool operator() (HitT const & a, HitT const & b) const {
	return ( a.time() < b.time() );
      }

  };


  // utility functor to sort hits by id and time
  template <typename HitT> class lessByIdAndTime {

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
  template <typename HitT> class lessByCIdAndTime {

     public:


       explicit lessByCIdAndTime(Calorimeter const & cal): _cal(0) {}

       bool operator() (HitT const & a, HitT const & b) const {
	  return ( _cal.crystalByRO(a.id()) < _cal.crystalByRO(b.id()) ||
        	  (_cal.crystalByRO(a.id()) == _cal.crystalByRO(b.id()) &&
                   a.time() < b.time()
                   )
        	  );
       }


     private:

       Calorimeter const & _cal;

  };


  // utility functor to sort hits by crystal id and time
  template <typename HitT> class lessByCIdAndTimeByPointer {

     public:

       explicit lessByCIdAndTimeByPointer(Calorimeter const * cal): _cal(cal) {}

       bool operator() (HitT const * a, HitT const * b) const {

	   return (  _cal->crystalByRO(a->id()) <  _cal->crystalByRO(b->id()) ||
        	   ( _cal->crystalByRO(a->id()) == _cal->crystalByRO(b->id()) &&
                     a->time() < b->time()
        	   )
        	  );		     
       }

     private:

       Calorimeter const * _cal;

  };

}

#endif /* Mu2eUtilities_sort_functors_hh */
