//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkErrCode.hh,v 1.13 2004/08/06 06:31:40 bartoldu Exp $
//
// Description:
//     Encapsulate error/success status of tracking operations.
//       Either failure() or success() will be non-zero, but not both.
//       Failure => no valid answer available.
//       Success => a valid answer has been
//       provided, even if it wasn't exactly what you asked for.  The
//       value of failure() or success() distinguishes different
//       failure/success modes.  A string describing the success/failure
//       mode can also be provided, and printed by the user.
//
//       Note that if this string is provided by the called function,
//       it _must_ be a pointer to a statically stored string (which includes
//       string literals).  E.g.
//           TrkErrCode err;
//           err.setFailure(10,"Forgot to tie my shoelaces.");
//           return err;
//       is valid.
//
//     Several codes have predefined meanings and strings; strings
//       supplied for them will be ignored.  Strings for codes >= 10
//       can be supplied by users.  Predefined:
//     failure = 1 -- "Arithmetic error."
//             = 2 -- "Failed to converge."
//             = 3 -- "Failed because parallel."
//             = 4-9 -- reserved until I think of some more standard codes
//     success = 1 -- "Normal completion."
//             = 2 -- "Didn't converge."
//             = 3 -- "Parallel"
//             = 4-9 -- reserved
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Authors: (Steve Schaffner) -- initial implementation stolen from A. Snyder
//
//------------------------------------------------------------------------
#ifndef TRKERRCODE_HH
#define TRKERRCODE_HH

#include <iosfwd>
#include <string>

// Class interface //
class TrkErrCode {
public:

  enum TrkSuccess {fail,succeed};

  TrkErrCode(TrkSuccess=succeed, int code=1, const char* str=0);
  ~TrkErrCode();

  // Copy constructors
  TrkErrCode(const TrkErrCode&);
  TrkErrCode& operator=(const TrkErrCode&);

  // access
  int failure() const                            {return _failure;}
  int success() const                            {return _success;}
  const std::string& message() const
  {
    return (_string != 0) ? *_string : _nullStr;
  }
  void print(std::ostream& ostr) const;

  // set
  void setMessage(const char* str = 0) {
    if (_string != 0) delete _string;
    if (str != 0) {
      _string= new std::string(str);
    }
    else {
      _string = 0;
    }
  }
  void setFailure(int i, const char* str = 0)
  {
    setMessage(str);
    _failure=(i==0?1:i); _success=0;
  }
  void setSuccess(int i, const char* str = 0) {
    setMessage(str);
    _success=(i==0?1:i); _failure=0;
  }


private:
  // Data
  int _failure;
  int _success;
  std::string*  _string;
  static std::string _nullStr;
};

std::ostream& operator<<(std::ostream& os, const TrkErrCode& trkerr);

#endif


