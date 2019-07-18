// The TrackCuts::accepted() method applies cuts on reconstructed
// tracks and fills histograms of reconstructed quantities.  No MC
// dependencies.
//
// Andrei Gaponenko, 2014

#ifndef Mu2eUtilities_inc_TrackCuts_hh
#define Mu2eUtilities_inc_TrackCuts_hh

#include <string>

#include "art_root_io/TFileDirectory.h"

#include "Mu2eUtilities/inc/HistTrackSum.hh"

namespace fhicl { class ParameterSet; }
class TH1;

namespace mu2e {
  class TrackSummary;

  class TrackCuts {
  public:

    bool accepted(const TrackSummary& trk, double weight);

    TrackCuts(const fhicl::ParameterSet& pset, art::TFileDirectory& topdir, const std::string& subdir);

  private:
    int cutnactive_;
    double cutfitcon_;
    double cutmomerr_;
    double cutt0err_;
    double cutd0min_;
    double cutd0max_;
    double cutdoutmin_;
    double cutdoutmax_;
    double cutt0min_;
    double cutt0max_;
    double cuttandipmin_;
    double cuttandipmax_;
    double cutmommin_;
    double cutmommax_;

    // A kludge to make sure our ROOT directory is created before
    // the nested directories for the HistTrackSum objects
    static art::TFileDirectory makeMyTopDir(art::TFileDirectory& parent, const std::string& subdir);
    art::TFileDirectory mytopdir_;

    HistTrackSum hall_;
    HistTrackSum hfinal_;
    TH1 *nactive_, *fitcon_, *momerr_, *t0err_, *t0_, *d0_, *dout_, *tanDip_, *momentum_;
  };
}

#endif/*Mu2eUtilities_inc_TrackCuts_hh*/
