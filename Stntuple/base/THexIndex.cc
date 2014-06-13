//

#include "base/THexIndex.hh"


THexIndex THexIndex::fgPos[6] = {
  THexIndex( 1, 0),  // up right
  THexIndex( 0, 1),  // down right
  THexIndex(-1, 1),  // down
  THexIndex(-1, 0),  // down left
  THexIndex( 0,-1),  // up left
  THexIndex( 1,-1)   // up 
};
