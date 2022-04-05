
//  Rigid body alignment parameters for tracker
//
#ifndef DbTables_TrkAlignParams_hh
#define DbTables_TrkAlignParams_hh

#include "Offline/GeneralUtilities/inc/HepTransform.hh"
#include "Offline/DataProducts/inc/StrawId.hh"
#include "Offline/DataProducts/inc/StrawIdMask.hh"
#include "cetlib_except/exception.h"

namespace mu2e {
  class TrkAlignParams {
    public:
      TrkAlignParams(int index, StrawId const& id, float dx, float dy, float dz, float rx, float ry, float rz):
        _index(index), _id(id), _dx(dx),_dy(dy),_dz(dz),
        _rx(rx),_ry(ry),_rz(rz),_transform(dx,dy,dz,rx,ry,rz) {}
        // the following should really use float, but the alignment classes want double FIXME!
      TrkAlignParams(StrawId const& id, StrawIdMask::Level level, double dx, double dy, double dz, double rx, double ry, double rz):
        _id(id), _dx(dx),_dy(dy),_dz(dz),
        _rx(rx),_ry(ry),_rz(rz),_transform(dx,dy,dz,rx,ry,rz) {
          switch(level) {
            case StrawIdMask::tracker:
              _index = 0;
              break;
            case StrawIdMask::plane:
              _index = id.plane();
              break;
            case StrawIdMask::uniquepanel:
              _index = id.uniquePanel();
              break;
            default:
              throw cet::exception("Geom") << "Illegal level: " << level << std::endl;
          }
        }

      int index() const { return _index; }
      StrawId const& id() const { return _id; }
      float dx() const {return _dx;}
      float dy() const {return _dy;}
      float dz() const {return _dz;}
      float rx() const {return _rx;}
      float ry() const {return _ry;}
      float rz() const {return _rz;}
      HepTransform const& transform() const {return _transform;}
    private:
      int _index;// index is redundant with StrawId and the table row count but is required by the Db interface
      StrawId _id;
      float _dx;
      float _dy;
      float _dz;
      float _rx;
      float _ry;
      float _rz;
      HepTransform _transform;
  };
}
#endif

