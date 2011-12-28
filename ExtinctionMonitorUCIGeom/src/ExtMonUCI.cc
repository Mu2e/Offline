// $Id: ExtMonUCI.cc,v 1.3 2011/12/28 00:25:05 youzy Exp $
// $Author: youzy $
// $Date: 2011/12/28 00:25:05 $

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
    const ExtMonMag* ExtMon::mag(unsigned int iMag) const
    {
      if (iMag < _mags.size()) return &_mags[iMag];
      else
      {
        std::cout << "ExtMonUCI::mag " << iMag << " >= " << "mags size " << _mags.size() << std::endl;
        return 0;
      }
    }

    //================================================================
    const ExtMonTof* ExtMon::tof(unsigned int iTofStation, unsigned int iTofSegment) const
    {
      unsigned int iTof = iTofStation * _nTofSegments + iTofSegment;
      if (iTof < _tofs.size()) return &_tofs[iTof];
      else
      {
        std::cout << "ExtMonUCI::tof " << iTof << " >= " << "tofs size " << _tofs.size() << std::endl;
        return 0;
      }
    }

    //================================================================
    const ExtMonTof* ExtMon::tof(unsigned int iTof) const
    {
      if (iTof < _tofs.size()) return &_tofs[iTof];
      else
      {
        std::cout << "ExtMonUCI::tof " << iTof << " >= " << "tofs size " << _tofs.size() << std::endl;
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
