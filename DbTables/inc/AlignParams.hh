//
//  Rigid body alignment parameters
//
#ifndef DbTables_AlignParams_hh
#define DbTables_AlignParams_hh

#include "GeneralUtilities/inc/HepTransform.hh"

namespace mu2e {
  class AlignParams {
    public:
      AlignParams(int index, float dx, float dy, float dz, 
	  float rx, float ry, float rz):
	_index(index),_dx(dx),_dy(dy),_dz(dz),
	_rx(rx),_ry(ry),_rz(rz),_transform(dx,dy,dz,rx,ry,rz) {}
      int index() const { return _index; }
      float dx() const {return _dx;}
      float dy() const {return _dy;}
      float dz() const {return _dz;}
      float rx() const {return _rx;}
      float ry() const {return _ry;}
      float rz() const {return _rz;}
      HepTransform const& transform() const {return _transform;}
    private:
      int _index;
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

