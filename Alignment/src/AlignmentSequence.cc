// AlignmentSequence.cc
// This is the definitions file for the AlignmentSequence class.
// Notes about the function of the class can be found in 
// Alignment/inc/AlignmentSequence.hh
// Original Author:  David Norvil Brown (UofLouisville), Oct 2017

#include "Alignment/inc/AlignmentSequence.hh"

namespace mu2e {

  AlignmentSequence::AlignmentSequence(const AlignmentSequence& rhs ) {
    _lastValidity = rhs.lastValidity();
    unsigned int count = rhs.size();
    for ( unsigned int iCopy = 0; iCopy < count; iCopy++ ) {
      std::pair<IoV,AlignmentObj> anAlign(rhs.at(iCopy).first,rhs.at(iCopy).second);
      _sequence.push_back(anAlign);
    }

  } // end of copy ctor


  void AlignmentSequence::addPair ( IoV & interval, AlignmentObj & alignObj ) {

    std::pair<IoV,AlignmentObj> aPr(interval,alignObj);
    _sequence.push_back(aPr);

  } // end of def of addPair

  AlignmentObj AlignmentSequence::getAlignment( const unsigned long& aTime ) {
    if ( mu2e::IoV::isValid(aTime) ) {
      bool found = false;
      if ( _sequence[_lastValidity].first.contains(aTime)) {
	found = true;
	return _sequence[_lastValidity].second;
      }
      int checkHi = (int)_lastValidity + 1; 
      int checkLo = (int)_lastValidity - 1;
      while ( !found && ( checkHi < (int)_sequence.size() || checkLo > -1 ) ) {
	if ( checkHi < (int)_sequence.size() ) {
	  if ( _sequence[checkHi].first.contains(aTime)) {
	    _lastValidity = checkHi;
	    return _sequence[checkHi].second;
	  }
	}
	if ( checkLo > -1 ) {
	  if ( _sequence[checkLo].first.contains(aTime)) {
	    _lastValidity = checkLo;
	    return _sequence[checkLo].second;
	  }
	}
	checkHi++;
	checkLo--;
      } // end of while, searching for appropriate IoV
    }
    // If it gets here, there is a problem.  Someone asked for an invalid time
    // Right now, just return a null AlignmentObj.  Will want to think about
    // right thing to do!
    AlignmentObj tmpNullAO;
    return tmpNullAO;
  }

  std::ostream& operator<<(std::ostream& os, const AlignmentSequence& rhs) {
    os << "AlignmentSequence:  {\n";
    for ( unsigned int i = 0; i < rhs.size(); i++ ) {
      os << rhs.at(i).first << ":  " << rhs.at(i).second << "\n";
    } // end of loop
    os << "\t }";
    return os;
  } // end of output operator

}  // end of namespace mu2e
