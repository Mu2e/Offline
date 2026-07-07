//
// Read the table with data about DIO spectrum, and
// merge the spectrum with the analytic expression
// in the endpoint region taken from Czarnecki spectrum
// Czarneckki et al 10.1103/PhysRevD.84.013006
//
//

#include <array>
#include <string>
#include <utility>

#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"
// Mu2e includes
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"
// Framework includes
#include "cetlib/pow.h"

#include "Offline/Mu2eUtilities/inc/CzarneckiSpectrum.hh"
#include "Offline/Mu2eUtilities/inc/Table.hh"

namespace mu2e {

  CzarneckiSpectrum::CzarneckiSpectrum() :
    _table ( loadTable<2>( ConfigFileLookupPolicy()( "Offline/EventGenerator/data/czarnecki_"+
                                                     GlobalConstantsHandle<PhysicsParams>()->getStoppingTargetMaterial()+".tbl" ) ) )
  {
    _halfBinWidth = (_table.getRow(1).first - _table.getRow(0).first)/2.;
  }

  double CzarneckiSpectrum::getWeight(double E) const {

    // FIXME remove bin centering
    // assume E is bin center, but table is listing bin left edge
    double E_leftedge = E - _halfBinWidth;

    const unsigned iRow = _table.getLowerBoundRowIndex( E_leftedge );

    double weight       = _table( iRow ) ;
    if ( iRow == _table.getNrows()-1 || iRow == 0 ) return weight;

    auto const & row        = _table.getRow( iRow   );
    auto const & row_before = _table.getRow( iRow-1 );
    auto const & row_after  = _table.getRow( iRow+1 );

    weight = interpolate( E_leftedge, row_after, row, row_before );
    if ( weight < 0 ) weight = interpolateE5 ( E_leftedge, row );

    return weight;
  }


  double CzarneckiSpectrum::interpolate (const double E,
                                         const TableRow<2>& row_after,
                                         const TableRow<2>& row,
                                         const TableRow<2>& row_before ) const {

    const double e1(  row_after.first );  const double p1(  row_after.second.at(0) );
    const double e2(        row.first );  const double p2(        row.second.at(0) );
    const double e3( row_before.first );  const double p3( row_before.second.at(0) );

    const double discr = e1*e1*e2 + e1*e3*e3 + e2*e2*e3 - e3*e3*e2 - e1*e1*e3 - e1*e2*e2;

    const double A = (p1*e2 + p3*e1 + p2*e3 - p3*e2 - p1*e3 - p2*e1) / discr;

    const double B = (e1*e1*p2 + e3*e3*p1 + e2*e2*p3 - e3*e3*p2 - e1*e1*p3 - e2*e2*p1) / discr;

    const double C = (e1*e1*e2*p3 + e3*e3*e1*p2 + e2*e2*e3*p1 -
                      e3*e3*e2*p1 - e1*e1*e3*p2 - e2*e2*e1*p3) / discr;

    return (A*E*E + B*E + C);

  }

  double CzarneckiSpectrum::interpolateE5(double E, TableRow<2> val ) const {

    GlobalConstantsHandle<PhysicsParams> phy;

    const double energy = val.first;
    const double weight = val.second.at(0);

    const double b = weight/cet::pow<5>( phy->getMuonEnergy()-energy-cet::square(energy)/(2*phy->getAtomicMass()));

    if ( phy->getEndpointEnergy() - E < 0 ) return 0.;
    else return b*cet::pow<5>( phy->getMuonEnergy() - E - cet::square(E)/(2*phy->getAtomicMass()) );
  }

}

