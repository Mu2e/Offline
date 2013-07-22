#ifndef Mu2eUtilities_CzarneckiSpectrum_hh
#define Mu2eUtilities_CzarneckiSpectrum_hh
//
// Read Czarnecki DIO spectrum from a table and merge it
// with the spectrum coming from the endopoint region formula

// $Id: CzarneckiSpectrum.hh,v 1.6 2013/07/22 18:57:42 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2013/07/22 18:57:42 $
//
// Original Author: Kyle Knoepfel
//

// Mu2e includes
#include "Mu2eUtilities/inc/Table.hh"

// C++ includes
#include <utility>

namespace mu2e {

  class CzarneckiSpectrum {

  public:
    
    CzarneckiSpectrum();

    double getWeight(double E);

  private:

    const Table<2> _table;

    double interpolate(const double E,
                       const TableRow<2>& row_after, 
                       const TableRow<2>& row,
                       const TableRow<2>& row_before);

    double interpolateE5(double E, TableRow<2> value);

  };

} // end of namespace mu2e

#endif /* Mu2eUtilities_CzarneckiSpectrum_hh */

