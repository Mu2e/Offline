#ifndef GeneralUtilties_PathnameWithNextVersion_hh
#define GeneralUtilties_PathnameWithNextVersion_hh
//
// Utility to help manage pathnames of the form:
//
//  (optional path prefix) (filename)  (.) (version number)
//
// where the the version number is string that contains
// only demical digits.  The filename may or may not
// contain an extension ( .txt, .log etc ); if present it
// is considered part of the filename.  The optional path
// prefix may be relative, abolute or missing.
//
// For purposes of this documentation we define two ideas:
//    1) Versioned pathname.  This is the structure given above
//    2) Base pathname.  This contains only the elements:
//         (optional path prefix) (filename)
//       Both the dot and the version number are not included.
//
// The contructor takes an argument that is a base pathname.
//
// The constructor will look for any existing files that
// match the format of the versioned pathname; it will
// return a pathname with a version number that is one
// greater than the largest existing version number.
//
// If there are no existing files that match the format
// of the versioned pathname, the code will return a
// pathname with a version of 1.
//
// Limitations:
//
// 1) This code does not attempt to deal with race conditions
//    that may occur if two jobs are both trying to write files
//    to the same directory.
//
// Notes:
//
//  1) The only reason that I made this a class, not a free function,
//     was to provide both the output pathname as a string and the
//     version number as an int. If we decide that the version number
//     accessor is not needed this is better as a free function.
//

#include <string>

class PathnameWithNextVersion{
public:

  explicit PathnameWithNextVersion ( std::string const& inputPathname );

  std::string const& pathname() const { return pathname_; }
  size_t              version() const { return version_;  }

private:

  // Output path name, with the next available version number included.
  std::string pathname_;

  // The next available version number.
  std::size_t version_ = 1;

};

#endif
