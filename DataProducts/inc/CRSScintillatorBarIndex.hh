#ifndef DataProducts_CRSScintillatorBarIndex_hh
#define DataProducts_CRSScintillatorBarIndex_hh

//
// Dense integer identifier of one CRSScintillatorBar.
// Has values 0...(N-1), where N is the number
// of CRSScintillatorBars in the system.

//
// $Id: CRSScintillatorBarIndex.hh,v 1.1 2012/03/29 22:59:13 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/03/29 22:59:13 $
//
// Original author KLG; based on Rob Kutschke StrawIndex
//

#include <ostream>

namespace mu2e {


  class CRSScintillatorBarIndex{

  public:

    // No default c'tor by design.

    // No automatic conversion of int to CRSScintillatorBarIndex.
    explicit CRSScintillatorBarIndex(int idx):
      _idx(idx){
    }

    // Compiler generated versions are OK for:
    // copy c'tor, destructor, operator=

    // Return the value as an int or as an unsigned in
    // Do not want automatic conversion to an int.
    int      asInt() const { return _idx;}
    unsigned asUint() const { return static_cast<unsigned>(_idx);}

    bool operator==( CRSScintillatorBarIndex const & rhs) const{
      return (_idx == rhs._idx);
    }

    bool operator<( CRSScintillatorBarIndex const & rhs) const{
      return ( _idx < rhs._idx);
    }
    bool operator>( CRSScintillatorBarIndex const & rhs) const{
      return ( _idx > rhs._idx);
    }
  private:

    int _idx;
  };

  inline std::ostream& operator<<( std::ostream& ost,
                                   CRSScintillatorBarIndex const & idx){
    ost << idx.asInt();
    return ost;
  }

  inline bool operator!=( CRSScintillatorBarIndex const & lhs,
                          CRSScintillatorBarIndex const & rhs) {
    return !( lhs == rhs);
  }

}
#endif /* DataProducts_CRSScintillatorBarIndex_hh */
