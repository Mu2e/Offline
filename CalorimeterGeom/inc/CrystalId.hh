#ifndef CRYSTALID_HH
#define CRYSTALID_HH
//
// Identifier of one crystal in a calorimeter.
//

//
// $Id: CrystalId.hh,v 1.10 2010/08/10 19:06:58 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/08/10 19:06:58 $
//
// Original author R. Bernstein and Rob Kutschke
//
#include <ostream>
#include "CalorimeterGeom/inc/RSliceId.hh"

namespace mu2e { 
  namespace calorimeter{
    struct CrystalId{

    public:

      
      /*      CrystalId():
              _rslid(RSliceId()),
              _n(-1){
              }

              CrystalId( RSliceId rslice,
              int n
              ):
              _rslid(rslice),
              _n(n){
              }

              CrystalId( ZSliceId zsliceid,
              int rslice,
              int n
              ):
              _rslid(zsliceid,rslice),
              _n(n){
              }


      */
      CrystalId():
        _rslid(RSliceId()){}
      
      
      CrystalId( RSliceId rslice):
        _rslid(rslice){}
  
      CrystalId( ZSliceId zsliceid,
                 int rslice
                 ):
        _rslid(zsliceid,rslice){}

      CrystalId( VaneId vane,
                 int zslice,
                 int rslice
                 ):
        _rslid(RSliceId(vane,zslice,rslice)) {}

      ~CrystalId  (){
      }

      // Compiler generated copy c'tor and assignment
      // operators should be should be OK.

      const VaneId& getVaneId() const {
        return _rslid._zid._vid;
      }

      const ZSliceId& getZSliceId() const {
        return _rslid._zid;
      }

      const RSliceId& getRSliceId() const {
        return _rslid;
      }
  
      const int getVane() const{
        return _rslid._zid._vid;
      }

      const int getZSlice() const{
        return _rslid._zid._zslice;
      }

      const int getRSlice() const{
        return _rslid._rslice;
      }

      const int getCrystalIndex() const{
        return 0;
      }

      bool operator==( CrystalId const& rhs) const{
        return ( _rslid == rhs._rslid && _n == rhs._n );
      }

      bool operator!=( CrystalId const& rhs) const{
        return !( *this == rhs);
      }


      RSliceId _rslid;
      ZSliceId _zid;
      int _n;
  
    };

    inline std::ostream& operator<<(std::ostream& ost, 
                                    const CrystalId& c ){
      ost << "Crystal Id: ("
          << c.getRSliceId() 
          << " )";
      return ost;
    }

  } //namespace calorimeter
}
#endif
