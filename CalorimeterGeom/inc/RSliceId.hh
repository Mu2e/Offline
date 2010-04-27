#ifndef RSLICEID_HH
#define RSLICEID_HH

//
// Identifier of one rslice in a calorimeter.
//

//
// $Id: RSliceId.hh,v 1.1 2010/04/27 18:22:04 rhbob Exp $
// $Author: rhbob $
// $Date: 2010/04/27 18:22:04 $
//
// Original author R. Bernstein and Rob Kutschke
//

#include <ostream>

#include "CalorimeterGeom/inc/ZSliceId.hh"

namespace mu2e {
  namespace calorimeter{
    struct RSliceId{

    public:

      RSliceId():
	_zid(ZSliceId()),
	_rslice(-1){
      }

      RSliceId( ZSliceId zslice,
		int rslice
		):
	_zid(zslice),
	_rslice(rslice){
      }

      RSliceId( VaneId vane,
		int zslice,
		int rslice
		):
	_zid(ZSliceId(vane,zslice)),
	_rslice(rslice){
      }

      ~RSliceId  (){
      }

      const VaneId getVaneId () const{
	return _zid._vid;
      }
      const ZSliceId getZSliceId () const{
	return _zid;
      }

      const int getVane () const{
	return _zid._vid;
      }

      const int getZSlice () const{
	return _zid._zslice;
      }

      const int getRSlice() const{
	return _rslice;
      }

      bool operator==(RSliceId const& rhs) const{
	return ( _zid == rhs._zid && _rslice == rhs._rslice );
      }

      bool operator!=(RSliceId const& rhs) const{
	return !( *this == rhs);
      }

  
      // Compiler generated copy and assignment constructors
      // should be OK.

      ZSliceId _zid;
      int _rslice;
  
    };

    inline std::ostream& operator<<(std::ostream& ost, 
				    const RSliceId& rsl ){
      ost << rsl._zid << " " << rsl._rslice;
      return ost;
    }

  } // namespace calorimeter
} //namespace mu2e

#endif
