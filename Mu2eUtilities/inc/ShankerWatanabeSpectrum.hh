#ifndef Mu2eUtilities_ShankerWatanabeSpectrum_hh
#define Mu2eUtilities_ShankerWatanabeSpectrum_hh
//
// Read Watanabe data about DIO spectrum from a table and merge it
// with the spectrum coming from the Shanker formula

//
//

// Mu2e includes
#include "Mu2eUtilities/inc/Table.hh"

namespace mu2e {

  class ShankerWatanabeSpectrum {

  public:

    ShankerWatanabeSpectrum();
    ~ShankerWatanabeSpectrum(){}

    double getWeight(double E) const;

  private:

    double _wanaEndPoint, _wanaEndPointVal, _norm;

    // built-in energy tolerance (in MeV) - not used currently!
    const double _tolerance = 0.0049;

    const Table<2> _table;    

    double evaluateShanker (double E) const;
    double evaluateWatanabe(double E) const;

    double interpolate(const double E, 
                       const TableRow<2>& row_after,
                       const TableRow<2>& row,
                       const TableRow<2>& row_before ) const;

  };

} // end of namespace mu2e

#endif /* Mu2eUtilities_ShankerWatanabeSpectrum_hh */

