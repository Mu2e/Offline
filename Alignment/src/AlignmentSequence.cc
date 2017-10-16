// AlignmentSequence.cc
// This is the definitions file for the AlignmentSequence class.
// Notes about the function of the class can be found in 
// Alignment/inc/AlignmentSequence.hh
// Original Author:  David Norvil Brown (UofLouisville), Oct 2017

#include "Alignment/inc/AlignmentSequence.hh"

namespace mu2e {

  void AlignmentSequence::addPair ( IoV & interval, AlignmentObj & alignObj ) {

    std::pair<IoV,AlignmentObj> aPr(interval,alignObj);
    _sequence.push_back(aPr);

  } // end of def of addPair

  AlignmentObj AlignmentSequence::getAlignment( double aTime ) {
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
}  // end of namespace mu2e
