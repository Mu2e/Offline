/**
 * 
 * Remove leading and trailing whitespace from a string.
 * It modifies the input string.
 *
 * $Id: TrimInPlace.cc,v 1.1 2009/09/30 22:57:47 kutschke Exp $
 * $Author: kutschke $ 
 * $Date: 2009/09/30 22:57:47 $
 *
 * Original author Rob Kutschke
 * 
 */
#include <string>

#include "Mu2eUtilities/inc/TrimInPlace.hh"

using namespace std;

namespace mu2e {

// Remove leading and trailing white space.
void TrimInPlace( string& s){

  // Index of one past the character begin tested for non-white.
  int idx(s.size());

  // Find last non-white character.
  string::const_reverse_iterator b = s.rbegin();
  string::const_reverse_iterator e = s.rend();
  while (b!=e){
    char c = *b;
    if ( !isspace(c) ){
      int nn = s.size()-idx;
      s.erase(idx,nn);
      break;
    }
    --idx;
    ++b;
  }

  // The entire string is white space so erase it all.
  if ( idx == 0 ){
    s.erase(0,s.size());
    return;
  }
  
  // Find first non-white character.
  for ( string::size_type i=0;
	i<s.size(); ++i ){
    char c = s[i];
    if ( !isspace(c) ){
      s.erase(0,i);
      break;
    }
  }
}

}
