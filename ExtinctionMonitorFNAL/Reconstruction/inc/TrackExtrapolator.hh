// Propagate tracks in ExtMonFNAL
//
// Andrei Gaponenko, 2012

#ifndef ExtinctionMonitorFNAL_Reconstruction_TrackExtrapolator_hh
#define ExtinctionMonitorFNAL_Reconstruction_TrackExtrapolator_hh

namespace mu2e {

  class ExtMonFNALTrkParam;

  namespace ExtMonFNAL {

    class ExtMon;

    class TrackExtrapolator {
    public:
      // Ownership is not passed. The object must be kept available while the extrapolator is in use.
      explicit TrackExtrapolator(const ExtMon *detectorGeom);

      void extrapolateStraightLine(double newz0, ExtMonFNALTrkParam *inout) const;

      // The input point must be outside the magnet at positive z.
      // The call updates the argument to move the track to the magnet
      // entrance plane.  All coordinates are in the upstream stack
      // system.
      void extrapolateToMagnet(ExtMonFNALTrkParam *inout) const;

      // Input pars in the upstream stack system.  The point must be
      // on the magnet entrance plane.  Output pars are given in the
      // downstream stack system, the point is on the magnet exit
      // plane.  Returns true and updates the parameters on success,
      // returns false and does not alter the argument if the
      // extrapolation can't be done.  Note that not all possible
      // constraints are taken into account, therefore it is possible
      // that some "bad" tracks will be wrongly moved through the
      // magnet.  The code should do the right thing for anything
      // resembling a signal tracks.  The false return is also
      // definite.
      bool extrapolateThroughMagnet(ExtMonFNALTrkParam *inout) const;

      // Extrapolates to the given spectrometer plane.  The input z0
      // should be outside of the magnet, input pars in the UP stack
      // system for z0>0 or in the DN stack for z0<0.  On success
      // returns true and sets output pars to in the appropriate
      // system.  On failure returns false, *inout is messed up.
      bool extrapolateToPlane(unsigned plane, ExtMonFNALTrkParam *inout) const;

    private:
      const ExtMon *extmon_;
    };

  } // namespace ExtMonFNAL
} // namespace mu2e

#endif/*ExtinctionMonitorFNAL_Reconstruction_TrackExtrapolator_hh*/
