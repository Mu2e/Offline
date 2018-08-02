// AlignmentObj.cc
// This is the definitions file for the AlignmentObj class.
// This class holds a single alignment object for Mu2e geometry.
// David Norvil Brown, UofL, Oct 2017

#include "Alignment/inc/AlignmentObj.hh"

namespace mu2e {

  AlignmentObj::AlignmentObj( const AlignmentObj& rhs ) {
    _displacement = rhs.displacement();
    _rotation = rhs.rotation();
    _isValid = rhs.isValid();
  } // end of copy constructor

  std::ostream& operator<<(std::ostream& os, const AlignmentObj& rhs) {
    os << "AlignmentObj: translate: " << rhs.displacement() << ", rotate: " << rhs.rotation() << "Validity:  " << rhs.isValid() << ")\n";
    return os;
  } // end of outputter

} // end of namespace mu2e
