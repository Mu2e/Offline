#ifndef Mu2eUtilities_DIOBase_hh
#define Mu2eUtilities_DIOBase_hh

//
// Base class to allow generic access to all the classes that define
// DIO momentum spectrum.
//
// $Id: DIOBase.hh,v 1.2 2011/05/17 15:41:36 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:41:36 $
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
