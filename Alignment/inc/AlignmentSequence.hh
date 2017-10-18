#ifndef Alignment_AlignmentSequence_HH
#define Alignment_AlignmentSequence_HH
// AlignmentSequence.hh
// This is the header file for the AlignmentSequence class.
// The AlignmentSequence class holds sets of AlignmentObj objects
// and associated Intervals of Validity for a given geometric object.
// The AlignmentObj does not know about the AlignmentSequence that
// holds it, nor the interval of validity for which it is valid.
// Likewise, the AlignmentSequence does not know which geometric object
// it is for.
// In this implementation, the AlignmentSequence is effectively a 
// std::vector of std::pair, the first element of the pair being
// the Interval of Validity and the second element being the AlignmentObj.
// Original Author:  David Norvil Brown (UofLouisville), Oct 2017.

#include "Alignment/inc/AlignmentObj.hh"
#include "Alignment/inc/IoV.hh"
#include "Alignment/inc/ShapeDetail.hh"
#include <vector>
#include <utility>
#include <ostream>

namespace mu2e {
  class AlignmentSequence {
  public:
    AlignmentSequence() { _lastValidity = 0; }
    AlignmentSequence( const AlignmentSequence& rhs ); // copy c-tor

    void addPair( IoV & interval, AlignmentObj& alignObj );
    AlignmentObj getAlignment( double atTime );  // How do we want to handle invalid times?

    unsigned int lastValidity() const  { return _lastValidity; }
    unsigned int size() const  { return _sequence.size(); }

    std::pair<IoV,AlignmentObj> at(unsigned int& i ) const { return _sequence[i];}
    bool isValid() const { return ( _sequence.size() > 0 ); }

  private:
    std::vector<std::pair<IoV,AlignmentObj> >   _sequence;
    unsigned int                                _lastValidity;

  }; // end class def

  std::ostream& operator<<(std::ostream& os, const AlignmentSequence& rhs);
} // end namespace mu2e

#endif  //  Alignment_AlignmentSequence_HH
