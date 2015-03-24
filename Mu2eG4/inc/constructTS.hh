#ifndef Mu2eG4_constructTS_hh
#define Mu2eG4_constructTS_hh
//
// Free function to create  Transport Solenoid
//
// $Id: constructTS.hh,v 1.5 2013/06/28 19:26:33 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2013/06/28 19:26:33 $
//
// Original author KLG
//

namespace mu2e {

  class VolumeInfo;
  class SimpleConfig;
  class Beamline;

  void constructTS         ( VolumeInfo const& p, SimpleConfig const& c);

  void constructCryostat   ( VolumeInfo const& p, SimpleConfig const& c, Beamline const& bl);
  void constructCoils      ( VolumeInfo const& p, SimpleConfig const& c, Beamline const& bl);
  void constructCAs        ( VolumeInfo const& p, SimpleConfig const& c, Beamline const& bl);
  void constructCollimators( VolumeInfo const& p, SimpleConfig const& c, Beamline const& bl);
  void constructDegrader   ( VolumeInfo const& p, SimpleConfig const& c, Beamline const& bl);
  void constructPbarWindow ( VolumeInfo const& p, SimpleConfig const& c, Beamline const& bl);

}

#endif /* Mu2eG4_constructTS_hh */
