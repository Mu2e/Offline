#ifndef RecoDataProducts_ProtonBunchTime_hh
#define RecoDataProducts_ProtonBunchTime_hh
//
//  Data product representing the time of the peak of the proton bunch current at the front face of
//  the production target WRT the DAQ clock (EventWindowMarker)
//  of this event.  This can be estimated from early 'beam flash' detector signals.
//  Original author: Dave Brown (LBNL) 10/1/2020
//
namespace mu2e {
  struct ProtonBunchTime {
    float pbtime_;  // estimated mean time the proton bunch reaches the target
    float pbterr_; // estimated error on the mean
  };

}
#endif
