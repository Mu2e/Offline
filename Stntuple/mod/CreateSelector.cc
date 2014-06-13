//-----------------------------------------------------------------------------
//  Jul 26 2001 P.Murat  : STNTUPLE utility routines
//  Feb 17 2002 E.Thomson: add handling of views to StntupleGetProcessName
//-----------------------------------------------------------------------------
#include "Stntuple/mod/StntupleUtilities.hh"

//_____________________________________________________________________________
// void StntupleSetProcessName(StorableObject* Obj, const char* ProcessName) {
//   // this is clearly a hack, but it is needed for embedding - otherwise
//   // one would have to write a lot of code!

//   std::string* s = (std::string*) (((char*) Obj) + sizeof(TObject) + 
// 				   sizeof(Id)+sizeof(edm::RCPID));
//   *s = ProcessName;
// }

//_____________________________________________________________________________
void StntupleGetProcessName(const char* String, 
			    char*       ProcessName, 
			    char*       Description,
                            char*       CollType) 
{
  // parse string with form: PROCESS@DESCRIPTION;COLL_TYPE
  // by default COLL_TYPE is not defined, if present, we expect it to be
  // either COLL or VIEW (want to be able to use views as well as collections)
  // also allow old descriptions (PROCESS:DESCRIPTION;COLL_TYPE)
  
  const char* remainder;
//-----------------------------------------------------------------------------
// What is the process name? (PROD or L3 or USRP or not specified)
//-----------------------------------------------------------------------------
  const char* ptr;
  int         nb;
  ptr = strstr(String,"@");
  if (ptr) {
    nb              = ptr-String;
    strncpy(ProcessName,String,nb);
    ProcessName[nb] = 0;
    remainder       = ptr+1;
  } 
  else {
//-----------------------------------------------------------------------------
// backward compatibility issue - old branches may have ":" in the name
//-----------------------------------------------------------------------------
    ptr = strstr(String,":");
    if (ptr) {
      nb              = ptr-String;
      strncpy(ProcessName,String,nb);
      ProcessName[nb] = 0;
      remainder       = ptr+1;
    } 
    else {
      ProcessName[0] = 0;
      remainder      = String;
    }
  }
//-----------------------------------------------------------------------------
// What is the object type? (COLL or VIEW or not specified)
//-----------------------------------------------------------------------------
  ptr = strstr(remainder,";");
  if (ptr) {
    nb              = ptr-remainder;
    strncpy(Description,remainder,nb);
    Description[nb] = 0;
    if (CollType) strcpy(CollType,ptr+1);
  } 
  else {
    strcpy(Description,remainder);
    if (CollType) CollType[0] = 0;
  }
}

