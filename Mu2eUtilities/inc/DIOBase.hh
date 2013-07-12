#ifndef Mu2eUtilities_DIOBase_hh
#define Mu2eUtilities_DIOBase_hh

//
// Base class to allow generic access to all the classes that define
// DIO momentum spectrum.
//
// $Id: DIOBase.hh,v 1.4 2013/07/12 17:17:38 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2013/07/12 17:17:38 $
//
// Original author Kyle Knoepfel 
//                 

// C++ includes
#include <string>
#include <utility>
#include <vector>

namespace mu2e {

  class DIOBase {

  public:

    DIOBase() {
    }

    virtual ~DIOBase() {
    }

    virtual double getWeight(double E) = 0;

  protected:

    virtual void readTable() {}
    virtual void checkTable() const;

    typedef std::pair<double,double> SpectrumValue;

    std::vector<SpectrumValue> _table;

  };

} // end of namespace mu2e

#endif /* Mu2eUtilities_DIOBase_hh */
