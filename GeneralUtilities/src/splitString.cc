#include "Offline/GeneralUtilities/inc/splitString.hh"

#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>

using namespace boost;

namespace mu2e {

StringVec splitString(std::string const& input,
            std::string const& delimeter, std::string const& quote,
            std::string const& escape, bool qtrim, bool keepEmpty) {

  StringVec sv;
  escaped_list_separator<char> els(escape, delimeter, quote);
  tokenizer<escaped_list_separator<char>> tok(input, els);
  for (tokenizer<escaped_list_separator<char>>::iterator sitr = tok.begin();
             sitr != tok.end(); ++sitr) {
    std::string ss(*sitr);
    if (qtrim) trim(ss);
    if( (!ss.empty()) || keepEmpty) sv.emplace_back(ss);
  }

  return sv;
}

}
