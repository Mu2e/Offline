#ifndef VANE_HH
#define VANE_HH

//
// Hold information about one vane in the calorimter.
//

//
// $Id: Vane.hh,v 1.3 2010/05/18 20:29:14 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/05/18 20:29:14 $
//
// Original author R, Bernstein and Rob Kutschke
//

#include <vector>
//#include <iostream> 

#include "CalorimeterGeom/inc/VaneId.hh"
#include "CalorimeterGeom/inc/ZSlice.hh"
//using namespace std;
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
        //        std::cout << "zslice in getRSlice of Vane.hh = " << rslid.getZSlice()<< endl;

        // 5/7/10 this is a zslice object, and I am getting it's RSlice --> go to ZSlice next
        return _zslices.at(rslid.getZSlice()).getRSlice(rslid);
      }

      const Crystal& getCrystal ( const CrystalId& cid ) const{
        return _zslices.at(cid.getZSlice()).getCrystal(cid);
      }

      // Formatted string embedding the id of the zslice.
      std::string name( std::string const& base ) const;

      // On readback from persistency, recursively recompute mutable members.
      void fillPointers ( const Calorimeter& calorimeter ) const;


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
