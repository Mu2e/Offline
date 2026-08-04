#ifndef DbTables_CalAlignParams_hh
#define DbTables_CalAlignParams_hh
//
//  Rigid body alignment parameters for calorimeter
//
namespace mu2e {
  class CalAlignParams {
    public:
      CalAlignParams(int index, float dx, float dy, float dz,
                     float rx, float ry, float rz) :
        _index(index), _dx(dx), _dy(dy), _dz(dz),
        _rx(rx), _ry(ry), _rz(rz)
      {}

      unsigned index() const { return _index; }
      float dx()  const { return _dx; }
      float dy()  const { return _dy; }
      float dz()  const { return _dz; }
      float rx()  const { return _rx; }
      float ry()  const { return _ry; }
      float rz()  const { return _rz; }

    private:
      unsigned _index;
      float _dx;
      float _dy;
      float _dz;
      float _rx;
      float _ry;
      float _rz;
  };
}
#endif
