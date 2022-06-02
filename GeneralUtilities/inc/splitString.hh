#ifndef GeneralUtilities_splitString_hh
#define GeneralUtilities_splitString_hh

#include "Offline/GeneralUtilities/inc/StringVec.hh"

namespace mu2e {

  // split a string at a delimiter. strings between the quote chars are
  // not split. chars that are escaped are not parsed.
  // breaking at newlines is delimeter="\n" (note that multiline strings
  // often end in a newline so the last split word is empty)
  // allowing single and double quoted strings is quote="\"'"
  // unescaped quote characters are only for parsing, so are removed
  // qtrim=true to trim whitespace in result, keepEmpty to keep empty fields
  StringVec splitString(std::string const& input,
          std::string const& delimeter=",", std::string const& quote="\"",
          std::string const& escape="\\", bool qtrim=true, bool keepEmpty=true);

}
#endif
