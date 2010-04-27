#ifndef ZSLICEID_HH
#define ZSLICEID_HH

//
// Identifier for a zslice.
//

//
// $Id: ZSliceId.hh,v 1.1 2010/04/27 18:21:52 rhbob Exp $
// $Author: rhbob $
// $Date: 2010/04/27 18:21:52 $
//
// Original author R. Bernstein and Rob Kutschke
//

#include <ostream>
#include "CalorimeterGeom/inc/VaneId.hh"

namespace mu2e { 
  namespace calorimeter{
    struct ZSliceId{

    public:

      ZSliceId():
	_vid(-1),
	_zslice(-1){
      }
  
      ZSliceId( VaneId vane,
		int zslice
		):
	_vid(vane),
	_zslice(zslice){
      }
  
      ~ZSliceId  (){
      }
  
      // Compiler generated copy and assignment constructors
      // should be OK.

      const int getVaneId() const {
	return _vid;
      }

      const int getVane() const {
	return _vid;
      }

      const int getZSlice() const {
	return _zslice;
      }

      bool operator==(ZSliceId const& rhs) const{
	return ( _vid == rhs._vid && _zslice == rhs._zslice );
      }

      bool operator!=(ZSliceId const& rhs) const{
	return !( *this == rhs);
      }
  
      VaneId _vid;
      int _zslice;
  
    };

    inline std::ostream& operator<<(std::ostream& ost, 
				    const ZSliceId& zsl ){
      ost << zsl._vid << " " << zsl._zslice;
      return ost;
    }

  } // namespace calorimeter
}  //namespace mu2e

#endif
