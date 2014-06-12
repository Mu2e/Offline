// A set of histograms of reconstructed quantities, which can be
// filled from a TrackSummary object.  No MC dependencies.
//
// Andrei Gaponenko, 2014

#ifndef Mu2eUtilities_inc_HistTrackSum_hh
#define Mu2eUtilities_inc_HistTrackSum_hh

#include <string>

namespace art { class TFileDirectory; }
class TH1;

namespace mu2e {
  class TrackSummary;

  class HistTrackSum {
  public:
    void fill(const TrackSummary& trk, double weight);

    // Books histograms in the newly created subdirectory, with the
    // name of the subdir given by the last arg.  This is relative to
    // the root TFileDirectory for the current module.
    HistTrackSum(art::TFileDirectory topdir, const std::string& subdir);

    // noncopyable because of the TH*
    HistTrackSum(const HistTrackSum&) = delete;
    HistTrackSum& operator=(const HistTrackSum&) = delete;
  private:
    TH1 *nactive_, *fitcon_, *momerr_, *t0err_, *t0_, *d0_, *dout_, *tanDip_, *momentum_;
    TH1 *radius_, *wavelength_, *costh_;
  };
}

#endif/*Mu2eUtilities_inc_HistTrackSum_hh*/
