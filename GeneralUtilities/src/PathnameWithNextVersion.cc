//
// Utility to help manage pathnames of the form:
//
//  (optional path prefix) (filename)  (.) (version number)
//
// See the header file for more details.
//
// Notes:
// 1) Here is how to read the regex:
//      a) start of string
//      b) the base pathname plus a dot character
//      c) one or more digits
//      d) end of string
//
// 2) If the file is in the current working directory there is
//    a special case.
//      a) Set the parent directory to "./"
//      b) Strip the leading "./" from each directory entry.
//

#include "GeneralUtilities/inc/PathnameWithNextVersion.hh"

#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <boost/type_index.hpp>

PathnameWithNextVersion::PathnameWithNextVersion( std::string const& basePathname ):
  pathname_(){

  // Regex used to identify pathnames that match the pattern of the
  // base pathname plus a version number. See note 1).
  boost::regex re( "^(" + basePathname + ".)(\\d+?)$");

  // Create a boost::path to the parent directory.
  boost::filesystem::path base(basePathname);
  boost::filesystem::path parent(base.parent_path());

  // Special case of the parent directory being cwd; see note 2a).
  size_t offset = 0;
  if ( parent.string().empty() ){
    parent = boost::filesystem::path(".");
    offset = 2;
  }

  // Scan the parent directory; select all files, if any, that match the
  // pattern of the versioned pathname; of these, find the maximum version number.
  int maxVersion{0};
  for ( auto const& entry : boost::filesystem::directory_iterator(parent) ){

    // See note 2b).
    if ( boost::regex_match(entry.path().string().substr(offset), re) ) {

      // The extension contains the leading dot; remove it.
      std::string ext{entry.path().extension().string().substr(1)};

      // By construction of the regexp, ext contains only demical digits.
      maxVersion = std::max( maxVersion, std::stoi(ext) );

    }
  }

  // Form the output filename.
  if ( maxVersion > 0 ){
    version_ = maxVersion+1;
  }
  pathname_ = basePathname + "." + std::to_string(version_);

}
