#ifndef Sandbox_XBarService_hh
#define Sandbox_XBarService_hh

//
//
// A test service that prints tracer ouptut. Used to study order of calls.
// The XbarSerivce depends on the FooService.  The XBarService is just
// a copy of the BarService but it comes after Foo in an alphabetic sort.
//
// $Id: XBarService.hh,v 1.1 2012/07/24 20:00:28 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/07/24 20:00:28 $
//
// Contact person Rob Kutschke
//

// Framework include files
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Run.h"

namespace mu2e {

  class XBarService {
  public:
    XBarService(const fhicl::ParameterSet&, art::ActivityRegistry&);
    ~XBarService();

    void preBeginRun( art::Run const &run);

    void poke() const;

  private:

  };

}

#endif /* Sandbox_XBarService_hh */
