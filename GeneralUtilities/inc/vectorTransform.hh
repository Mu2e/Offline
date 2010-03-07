#ifndef GeneralUtilities_vectorTransform_hh
#define GeneralUtilities_vectorTransform_hh

//
// Copy a vector of type T to a vector of type U.
// This will work provided the following works:
//  T t;
//  U u = static_cast<U> t;
// 
//
// $Id: vectorTransform.hh,v 1.1 2010/03/07 21:59:33 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/03/07 21:59:33 $
//

#include <vector>

template<class T, class U>
void vectorTransform( const std::vector<T>& in,
		      std::vector<U>&       out){
  out.clear();
  out.reserve(in.size());
  for ( size_t i=0; i<in.size(); ++i ){
    out.push_back(static_cast<U> (in[i]) );
  }
  
}

#endif
