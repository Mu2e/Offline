
#ifndef CRYSTALID_HH
#define CRYSTALID_HH
//
// Identifier of one crystal in a calorimeter.
//

//
// $Id: CrystalId.hh,v 1.7 2010/05/18 20:29:10 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/05/18 20:29:10 $
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

      CrystalId( VaneId vane,
                 int section,
                 int rslice,
                 int n
                 ):
        _rslid(RSliceId(vane,section,rslice)),
        _n(n){
      }

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

      const int getCrystal() const{
        return _n;
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
          << c.getRSliceId() << " "
          << c._n
          << " )";
      return ost;
    }

  } //namespace calorimeter
}
#endif
