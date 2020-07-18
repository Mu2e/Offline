#ifndef GeneralUtilities_inc_to_string_hh
#define GeneralUtilities_inc_to_string_hh

//
// A collection of to_string functions for various types
// that are not part of the Mu2e code base.  For types within
// the Mu2e code base, put the function with the .hh and .cc files
// for that type.
//
// Rob Kutschke, 2020


#include "canvas/Persistency/Provenance/SubRunID.h"

#include <string>

namespace mu2e{

//  Create a string representing an art::SubrunID in the format
//     rrr<sep>sss
// where rrr is the run number,  sss is the subrun number and
// <sep> is a user supplied separator that defaults to '_'.
// If the optional arguments lr and lsr are non-zero the rrr
// and sss fields are zero padded to lengths lr and lsr.
// If lr and lsr are too short, they are ignored and the strings
// have their natural length.
//
// The resulting string can be used for naming files, ROOT directories
// etc that need an embedded SubRunID.
//
// Note that art::SubRunID has a stream insertion operator
// but it's out is designed to be read by people and has a very different format.

  std::string to_string( art::SubRunID const& id, int lr=0, int lsr=0, std::string const& separator = "_" );
}

#endif /* GeneralUtilities_inc_to_string_hh */
