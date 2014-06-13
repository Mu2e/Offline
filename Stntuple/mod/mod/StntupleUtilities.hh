//-----------------------------------------------------------------------------
//  Jan 11 2001 P.Murat: declarations for a set of STNTUPLE utility routines
//-----------------------------------------------------------------------------
#ifndef StntupleUtilities_hh
#define StntupleUtilities_hh

#include "Stntuple/mod/InitStntupleDataBlocks.hh"

void StntupleGetProcessName     (const char* string, 
				 char*       ProcessName, 
				 char*       Description,
				 char*       CollType = 0);

// void StntupleSetProcessName     (StorableObject* Obj, const char* ProcessName);

#endif
