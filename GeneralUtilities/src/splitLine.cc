//
// Split a string, s, into pieces, using the given delimiter.
//
// Original author Rob Kutschke
//

#include "GeneralUtilities/inc/splitLine.hh"

using namespace std;

namespace mu2e {

  void splitLine ( std::string const& s,
                   std::string const& delimiter,
                   std::vector<std::string>& parts){

    size_t len = delimiter.size();

    if ( len == 0 ){
      parts.push_back(s);
      return;
    }

    size_t i = 0;
    while ( true ){
      size_t i1 = s.find( delimiter, i );
      if ( i1 == string::npos ){
        size_t npos = (i==0) ? s.size() : s.size()-1;
        string tmp = s.substr(i,npos);
        if ( !tmp.empty() ) parts.push_back(tmp);
        break;
      }
      string tmp = s.substr(i,i1-i);
      parts.push_back(tmp);
      i = i1+len;
    }
  }


} // end namespace mu2e
