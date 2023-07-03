#include "Offline/CaloConditions/inc/CalCalib.hh"

#include <sstream>

using namespace std;

namespace mu2e {

  void CalCalib::print( ostream& out) const{
    out << "ADCMeV "
      <<  ( _calpar->ADC2MeV() ) << endl;
    out << "algID "
      <<  ( _calpar->algID() ) << endl;
    out << " time offset " <<  (_calpar->timeOffset()) << endl;
  }

} // namespace mu2e
