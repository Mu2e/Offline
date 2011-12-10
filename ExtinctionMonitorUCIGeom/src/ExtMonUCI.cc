#include "ExtinctionMonitorUCIGeom/inc/ExtMonUCI.hh"

namespace mu2e {

  namespace ExtMonUCI {

    //================================================================
    ExtMon::ExtMon(const std::vector<double>& envelopeParams, 
		   const std::vector<double>& envelopeOrigin) 
      
      : _envelopeParams(envelopeParams)
      , _envelopeOrigin(envelopeOrigin)
    {
    }

    //================================================================
    CLHEP::Hep3Vector ExtMon::origin() const
    {
      return CLHEP::Hep3Vector(_envelopeOrigin[0], _envelopeOrigin[1], _envelopeOrigin[2]);
    }

    //================================================================
    CLHEP::Hep3Vector ExtMon::originLocal() const
    {
      return origin() - _hallOriginInMu2e;
    }

    //================================================================
    const ExtMonCol* ExtMon::col(unsigned int iCol) const
    {
      if (iCol < _cols.size()) return &_cols[iCol];
      else 
      {
        std::cout << "ExtMonUCI::col " << iCol << " >= " << "cols size " << _cols.size() << std::endl;
        return 0;
      }
    }

    //================================================================
    CLHEP::Hep3Vector ExtMon::extMonToMu2ePoint( CLHEP::Hep3Vector const& v ) const
    {
      return v + origin();
    }

    //================================================================
    CLHEP::Hep3Vector ExtMon::mu2eToExtMonPoint( CLHEP::Hep3Vector const& v ) const
    {
      return v - origin();
    }

  }
}
