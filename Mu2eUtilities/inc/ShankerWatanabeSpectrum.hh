#ifndef Mu2eUtilities_ShankerWatanabeSpectrum_hh
#define Mu2eUtilities_ShankerWatanabeSpectrum_hh
//
// Read Watanabe data about DIO spectrum from a table and merge it
// with the spectrum coming from the Shanker formula

// $Id: ShankerWatanabeSpectrum.hh,v 1.2 2013/07/22 18:57:42 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2013/07/22 18:57:42 $
//
//

// Mu2e includes
#include "Mu2eUtilities/inc/Table.hh"

namespace mu2e {

  class ShankerWatanabeSpectrum {

  public:

    ShankerWatanabeSpectrum();
    ~ShankerWatanabeSpectrum(){}

    double getWeight(double E);

  private:

    double _wanaEndPoint, _wanaEndPointVal, _norm;

    // built-in energy tolerance (in MeV) - not used currently!
    const double _tolerance = 0.0049;

    const Table<2> _table;    

    double evaluateShanker (double E);
    double evaluateWatanabe(double E);

    double interpolate(const double E, 
                       const TableRow<2>& row_after,
                       const TableRow<2>& row,
                       const TableRow<2>& row_before );

  };

} // end of namespace mu2e

#endif /* Mu2eUtilities_ShankerWatanabeSpectrum_hh */

