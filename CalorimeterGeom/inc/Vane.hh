#ifndef VANE_HH
#define VANE_HH

//
// Hold information about one vane in the calorimter.
//

//
// $Id: Vane.hh,v 1.1 2010/04/27 18:22:30 rhbob Exp $
// $Author: rhbob $
// $Date: 2010/04/27 18:22:30 $
//
// Original author R, Bernstein and Rob Kutschke
//

#include <vector>

#include "CalorimeterGeom/inc/VaneId.hh"
#include "CalorimeterGeom/inc/ZSlice.hh"

namespace mu2e {
  namespace calorimeter{
    class Vane{

      friend class Calorimeter;
      friend class CalorimeterMaker;

    public:

      // A free function, returning void, that takes a const Vane& as an argument.
      typedef void (*VaneFunction)( const Vane& s);

      Vane():_id(-1){}
      Vane( const VaneId& id ):_id(id){}
      ~Vane(){}
 
      // Compiler generated copy and assignment constructors
      // should be OK.
  
      const VaneId Id() const { return _id;}

      const std::vector<ZSlice>& getZSlices () const{ 
	return _zslices;
      }

      const ZSlice& getZSlice ( int n) const { 
	return _zslices.at(n);
      }

      const ZSlice& getZSlice ( const ZSliceId& zid ) const{
	return _zslices.at(zid._zslice);
      }

      const RSlice& getRSlice ( const RSliceId& rslid ) const{
	return _zslices.at(rslid.getZSlice()).getRSlice(rslid);
      }

      const Crystal& getCrystal ( const CrystalId& cid ) const{
	return _zslices.at(cid.getZSlice()).getCrystal(cid);
      }

      // Formatted string embedding the id of the zslice.
      std::string name( std::string const& base ) const;


#ifndef __CINT__

      // Loop over all crystals and call F.
      // F can be a class with an operator() or a free function.
      template <class F>
      inline void Vane::forAllCrystals ( F& f) const{
	for ( std::vector<ZSlice>::const_iterator i=_zslices.begin(), e=_zslices.end();
	      i !=e; ++i){
	  i->forAllCrystals(f);
	}
      }

      // Loop over all crystals and call F.
      // F can be a class with an operator() or a free function.
      template <class F>
      inline void Vane::forAllRSlices ( F& f) const{
	for ( std::vector<ZSlice>::const_iterator i=_zslices.begin(), e=_zslices.end();
	      i !=e; ++i){
	  i->forAllRSlices(f);
	}
      }

      template <class F>
      inline void Vane::forAllZSlices ( F& f) const{
	for ( std::vector<ZSlice>::const_iterator i=_zslices.begin(), e=_zslices.end();
	      i !=e; ++i){
	  f(*i);
	}
      }

#endif

    protected:

      VaneId _id;
      std::vector<ZSlice> _zslices;

    };
  } // namespace calorimeter
} //namespace mu2e

#endif
