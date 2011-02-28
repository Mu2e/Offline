#ifndef DIOBASE_HH
#define DIOBASE_HH

//
// Base class to allow generic access to all the classes that define
// DIO momentum spectrum.
//
// $Id: DIOBase.hh,v 1.1 2011/02/28 16:18:50 onoratog Exp $
// $Author: onoratog $
// $Date: 2011/02/28 16:18:50 $
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

#endif
