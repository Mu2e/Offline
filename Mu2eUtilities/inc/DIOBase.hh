#ifndef Mu2eUtilities_DIOBase_hh
#define Mu2eUtilities_DIOBase_hh

//
// Base class to allow generic access to all the classes that define
// DIO momentum spectrum.
//
// $Id: DIOBase.hh,v 1.3 2011/05/18 02:27:18 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:18 $
//
// Original author Gianni Onorato
//

namespace mu2e {

  class DIOBase {

  public:

    DIOBase() {
    }

    virtual ~DIOBase() {
    }

    virtual double fire() = 0;

  };

} // end of namespace mu2e

#endif /* Mu2eUtilities_DIOBase_hh */
