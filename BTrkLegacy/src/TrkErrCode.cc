//--------------------------------------------------------------------------
//// File and Version Information:
////      $Id: TrkErrCode.cc,v 1.13 2004/08/06 06:31:40 bartoldu Exp $
////
//// Description:
////
////
//// Environment:
////      Software developed for the BaBar Detector at the SLAC B-Factory.
////
//// Author(s): Steve Schaffner
////
//// Revision History:
////  20000420  M. Kelsey -- Remove terminating endl in print().
////------------------------------------------------------------------------
#include <iostream>
#include <algorithm>
#include "Offline/BTrkLegacy/inc/TrkErrCode.hh"

using std::ostream;

std::string TrkErrCode::_nullStr("");

TrkErrCode::TrkErrCode(TrkSuccess succ, int code, const char* str)
  : _string(0)
{
  setMessage(str);
  if (succ) {
    _failure = 0; _success = code;
  } else {
    _success = 0; _failure = code;
  }
}


TrkErrCode::TrkErrCode(const TrkErrCode& theCode)
  : _failure(theCode._failure)
    , _success(theCode._success)
{
  if (theCode._string != 0) {
    _string = new std::string(*theCode._string);
  }
  else {
    _string = 0;
  }
}


TrkErrCode::~TrkErrCode() {
  if (_string != 0) {
    delete _string;
    _string = 0;
  }
}

  TrkErrCode&
TrkErrCode::operator =(const TrkErrCode& theCode)
{
  _failure = theCode._failure;
  _success = theCode._success;

  if (theCode._string != 0) {
    if (_string != 0) {
      *_string = *theCode._string;
    }
    else {
      _string = new std::string(*theCode._string);
    }
  }
  else {
    if (_string != 0) delete _string;
    _string = 0;
  }

  return *this;
}


void
TrkErrCode::print(ostream& ostr) const
{
  const char* pstatus = 0;
  int code;
  if (success()) {
    pstatus = "succeeded";
    code = success();
  } else {
    pstatus = "failed";
    code = failure();
  }
  std::string pstring;
  static const std::string failed[4] = { "Arithmetic error.",
    "Failed to converge.",
    "Failed because parallel.",
    "Undefined error." };

  static const std::string succeeded[4] = { "Normal completion.",
    "Didn't converge.",
    "Parallel.",
    "Undefined success state."};

  if (code > 0 && code < 10) {
    if (failure()) {
      pstring = failed[std::min(code-1,3) ];
    } else if (success()) {
      pstring = succeeded[std::min(code-1,3) ];
    }
  } else if (_string == 0 ) {
    pstring = "Unknown error.";
  } else {
    pstring = *_string;
  }

  ostr << "TrkErrCode: " << pstatus << ", code " << code
    << ".  Status: " << pstring.c_str();
}

  ostream&
operator<<(ostream& os, const TrkErrCode& trkerr)
{
  trkerr.print(os);
  return os;
}
