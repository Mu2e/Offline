#ifndef BOX_HH
#define BOX_HH

namespace mu2e {

  class Box {
  public:

    Box( double x, double y, double z) :
      _xhl(x), _yhl(y), _zhl(z){}

    double getXhalfLength() const { return _xhl; }
    double getYhalfLength() const { return _yhl; }
    double getZhalfLength() const { return _zhl; }

  private:

    double _xhl;
    double _yhl;
    double _zhl;
  };

}

#endif
