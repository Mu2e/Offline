///////////////////////////////////////////////////////////////////////////////
//
// fData[0] : cluster X 
// fData[1] : cluster Y 
// fData[2] : seedIndex of the position which returns the best helix
//            for the patterrecognition
// level
///////////////////////////////////////////////////////////////////////////////

#include "CalPatRec/inc/THackData.hh"

//-----------------------------------------------------------------------------
THackData::THackData(): TNamed("","") {
}

//-----------------------------------------------------------------------------
THackData::THackData(const char* Name, const char* Title): TNamed(Name,Title) {
}

//-----------------------------------------------------------------------------
THackData::~THackData() {
}
