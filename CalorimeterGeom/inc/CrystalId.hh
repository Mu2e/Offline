#ifndef CRYSTALID_HH
#define CRYSTALID_HH
//
// Identifier of one crystal in a calorimeter.
//

//
// $Id: CrystalId.hh,v 1.5 2010/04/27 18:21:12 rhbob Exp $
// $Author: rhbob $
// $Date: 2010/04/27 18:21:12 $
//
// Original author R. Bernstein and Rob Kutschke
//
#include <ostream>
#include "CalorimeterGeom/inc/RSliceId.hh"

namespace mu2e { 
  namespace calorimeter{
    struct CrystalId{

    public:

      CrystalId():
	_rid(RSliceId()),
	_n(-1){
      }
  
      CrystalId( RSliceId rslice,
		 int n
		 ):
	_rid(rslice),
	_n(n){
      }
  
      CrystalId( ZSliceId zsliceid,
		 int rslice,
		 int n
		 ):
	_rid(zsliceid,rslice),
	_n(n){
      }

      CrystalId( VaneId vane,
		 int section,
		 int rslice,
		 int n
		 ):
	_rid(RSliceId(vane,section,rslice)),
	_n(n){
      }

      ~CrystalId  (){
      }

      // Compiler generated copy c'tor and assignment
      // operators should be should be OK.

      const VaneId& getVaneId() const {
	return _rid._zid._vid;
      }

      const ZSliceId& getZSliceId() const {
	return _rid._zid;
      }

      const RSliceId& getRSliceId() const {
	return _rid;
      }
  
      const int getVane() const{
	return _rid._zid._vid;
      }

      const int getZSlice() const{
	return _rid._zid._zslice;
      }

      const int getRSlice() const{
	return _rid._rslice;
      }

      const int getCrystal() const{
	return _n;
      }

      bool operator==( CrystalId const& rhs) const{
	return ( _rid == rhs._rid && _n == rhs._n );
      }

      bool operator!=( CrystalId const& rhs) const{
	return !( *this == rhs);
      }


      RSliceId _rid;
      int _n;
  
    };

    inline std::ostream& operator<<(std::ostream& ost, 
				    const CrystalId& c ){
      ost << "Crystal Id: ("
	  << c.getRSliceId() << " "
	  << c._n
	  << " )";
      return ost;
    }

  } //namespace calorimeter
}
#endif
