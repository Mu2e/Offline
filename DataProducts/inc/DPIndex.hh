#ifndef DataProducts_DPIndex_hh
#define DataProducts_DPIndex_hh

//
// A persistable index into another data product.
//
//
// $Id: DPIndex.hh,v 1.2 2011/06/01 14:55:55 greenc Exp $
// $Author: greenc $
// $Date: 2011/06/01 14:55:55 $
//

// Framework includes.
#include "art/Persistency/Provenance/ProductID.h"

#include <ostream>

namespace mu2e {
  struct DPIndex{

    // The actual data for this struct.
    art::ProductID id;
    unsigned       index;

    DPIndex():
      id(),
      index(0){}

    DPIndex( art::ProductID const& id_, int index_):
      id(id_),
      index(index_){
    }

    // Compiler generated versions are OK for:
    // destructor, copy c'tor, assignment operator.

  };

  inline bool operator==(const DPIndex& lhs,
                         const DPIndex& rhs){
      return ( lhs.id == rhs.id && lhs.index == rhs.index );
  }

  inline bool operator!=(const DPIndex& lhs,
                         const DPIndex& rhs){
    return !(lhs==rhs);
  }

  // Sort first on ProductID and then on index.
  inline bool operator<(const DPIndex& lhs,
                        const DPIndex& rhs){
    return ( lhs.id < rhs.id ) ||
      ( lhs.id == rhs.id && lhs.index < rhs.index );
  }

  // ProductID does not define operators >, <=, >= so we would need
  // to fix that before defining those operators for this class.


  inline std::ostream& operator<<( std::ostream& ost,
                                   DPIndex const& dpi ){
    ost << "(" << dpi.id
        << ","  << dpi.index
        << ")";
    return ost;
  }


}

#endif /* DataProducts_DPIndex_hh */
