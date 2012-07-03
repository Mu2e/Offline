#ifndef RecoDataProducts_LocalHelixPayload_hh
#define RecoDataProducts_LocalHelixPayload_hh
//
// The persistent data associated with one locally valid helical element
// of a Kalman trajectory.
//
// $Id: LocalHelixPayload.hh,v 1.1 2012/07/03 03:27:24 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/07/03 03:27:24 $
//
// Contact person Rob Kutschke
//

// Mu2e

// C++ includes
#include <iosfwd>
#include <vector>

#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/SymMatrix.h"

namespace mu2e {

  class LocalHelixPayload{

  public:

    CLHEP::HepVector    const& params() const { return params_; }
    CLHEP::HepSymMatrix const& cov()    const { return cov_;    }

    double lowRange()    const { return lowRange_;    }
    double highRange()   const { return highRange_;   }
    double globalRange() const { return globalRange_; }

    LocalHelixPayload():
      params_(),
      cov_(),
      lowRange_(0),
      highRange_(0),
      globalRange_(0){
    }
    // Accept compiler written d'tor, copy c'tor and copy assignment.

    LocalHelixPayload( CLHEP::HepVector    const& aparams,
                       CLHEP::HepSymMatrix const& acov,
                       double                     alowRange,
                       double                     ahighRange,
                       double                     aglobalRange):
      params_(aparams),
      cov_(acov),
      lowRange_(alowRange),
      highRange_(ahighRange),
      globalRange_(aglobalRange){
    }

    // Print contents of the object.
    void print( std::ostream& ost, bool doEndl=true) const;

  private:

    CLHEP::HepVector    params_;
    CLHEP::HepSymMatrix cov_;

    double              lowRange_;
    double              highRange_;
    double              globalRange_;

  };

  inline std::ostream& operator<<( std::ostream& ost,
                                   LocalHelixPayload const& hit){
    hit.print(ost,false);
    return ost;
  }


} // namespace mu2e

#endif /* RecoDataProducts_LocalHelixPayload_hh */
