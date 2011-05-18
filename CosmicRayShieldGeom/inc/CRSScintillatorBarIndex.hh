#ifndef CosmicRayShieldGeom_CRSScintillatorBarIndex_hh
#define CosmicRayShieldGeom_CRSScintillatorBarIndex_hh

//
// Dense integer identifier of one CRSScintillatorBar.
// Has values 0...(N-1), where N is the number
// of CRSScintillatorBars in the system.

//
// $Id: CRSScintillatorBarIndex.hh,v 1.3 2011/05/18 02:27:15 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:15 $
//
// Original author KLG; based on Rob Kutschke StrawIndex
//

#include <ostream>

namespace mu2e {


  class CRSScintillatorBarIndex{

  public:

    // No default c'tor by design.

    // No automatic conversion of int to CRSScintillatorBarIndex.
    explicit CRSScintillatorBarIndex(int32_t idx):
      _idx(idx){
    }

    // Compiler generated versions are OK for:
    // copy c'tor, destructor, operator=

    // Return the value as an int or as an unsigned in
    // Do not want automatic conversion to an int.
    int32_t  asInt() const { return _idx;}
    uint32_t asUint() const { return static_cast<uint32_t>(_idx);}

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

    int32_t _idx;
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
#endif /* CosmicRayShieldGeom_CRSScintillatorBarIndex_hh */
