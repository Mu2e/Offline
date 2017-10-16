// AlignmentMap.cc
// This is the definitions file for the AlignmentMap class.
// This class holds a map from GeomHandles to corresponding AlignmentSequences
// for those GeomHandles that have them.
// David Norvil Brown, Oct 2017

#include "Alignment/inc/AlignmentMap.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include <vector>
#include <string>
#include <iostream>

namespace mu2e {

  void AlignmentMap::make( const SimpleConfig& _config ) {

    // Get all the variable names found
    std::vector<std::string> theNames; 
    _config.getNames( theNames );
    for ( unsigned int iName = 0; iName < theNames.size(); iName++ ) {
      std::cout << theNames[iName] << std::endl;
    }

  } // end of make function

  AlignmentSequence AlignmentMap::find(std::string& handle) const {
    std::cout << "The handle is: " << handle << std::endl;
    AlignmentSequence tmpSeq;
    return tmpSeq;
  } // end of alignmentS

  void AlignmentMap::addAlignment(std::string& handle, AlignmentSequence& as) {
    _map.emplace(handle,as);
  }// end of addAlignment function def

} // end of namespace mu2e
