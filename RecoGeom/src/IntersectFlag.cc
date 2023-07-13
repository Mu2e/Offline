#include "Offline/RecoGeom/inc/IntersectFlag.hh"
#include <stdexcept>
#include <iostream>
#include <stdio.h>

std::ostream& operator <<(std::ostream& ost, mu2e::RecoGeom::IntersectFlag const& iflag) {
  if(iflag.onsurface_){
    ost << "On surface ";
    if(iflag.inbounds_){
      ost << " and in bounds ";
    }
  }else {
    ost << "No Intersection";
  }
  return ost;
}
