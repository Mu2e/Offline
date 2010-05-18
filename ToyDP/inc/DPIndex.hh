#ifndef DPIndex_HH
#define DPIndex_HH

//
// A persistable index into another data product.
//
//
// $Id: DPIndex.hh,v 1.5 2010/05/18 22:06:46 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/05/18 22:06:46 $
//

// Framework includes.
#include "DataFormats/Provenance/interface/ProductID.h"

namespace mu2e {
  struct DPIndex{

    // The actual data for this struct.
    edm::ProductID id;
    uint32_t       index;

    DPIndex():
      id(),
      index(0){}

    DPIndex( edm::ProductID const& id_, int index_):
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

#endif
