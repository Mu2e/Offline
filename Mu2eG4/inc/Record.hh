#ifndef RECORD_HH
#define RECORD_HH

#include "DiskRecord.hh"

struct Record{

  Record():
    x(),
    y(),
    z(),
    bx(),
    by(),
    bz(){
  }

  Record(DiskRecord const& r):
    x(r.x),
    y(r.y),
    z(r.z),
    bx(r.bx),
    by(r.by),
    bz(r.bz){
  }

  Record(double& x, double& y, double& z, double bx, double by, double bz):
    x(x),
    y(y),
    z(z),
    bx(bx),
    by(by),
    bz(bz){
  }

  ~Record(){}

  float x;
  float y;
  float z;
  float bx;
  float by;
  float bz;

};

#endif
