#include <exception>
#include <memory>
#include <vector>

#include "CLHEP/Vector/ThreeVector.h"

#include "RecoDataProducts/inc/TrackSummary.hh"
#include "TH1.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/exception.h"

#include "Mu2eUtilities/inc/TrackCuts.hh"

namespace mu2e {

  art::TFileDirectory TrackCuts::makeMyTopDir(art::TFileDirectory& parent, const std::string& subdir) {
    return subdir.empty() ? parent : parent.mkdir(subdir.c_str());
  }

  TrackCuts::TrackCuts(const fhicl::ParameterSet& pset,
                       art::TFileDirectory& topdir,
                       const std::string& subdir)
    : cutnactive_(pset.get<int>("nactive"))
    , cutfitcon_(pset.get<double>("fitcon"))
    , cutmomerr_(pset.get<double>("momerr"))
    , cutt0err_(pset.get<double>("t0err"))
    , cutd0min_(pset.get<double>("d0min"))
    , cutd0max_(pset.get<double>("d0max"))
    , cutdoutmin_(pset.get<double>("doutmin"))
    , cutdoutmax_(pset.get<double>("doutmax"))
    , cutt0min_(pset.get<double>("t0min"))
    , cutt0max_(pset.get<double>("t0max"))
    , cuttandipmin_(pset.get<double>("tandipmin"))
    , cuttandipmax_(pset.get<double>("tandipmax"))
    , cutmommin_(pset.get<double>("mommin"))
    , cutmommax_(pset.get<double>("mommax"))

    , mytopdir_(makeMyTopDir(topdir,subdir))
    , hall_(mytopdir_, "all")
    , hfinal_(mytopdir_, "final")
  {
    art::TFileDirectory dir = mytopdir_.mkdir("nm1", "Distributions after N-1 cut");
    TH1::SetDefaultSumw2();
    nactive_ = dir.make<TH1D>("nactive", "nactive", 150, -0.5, 149.5);
    fitcon_ = dir.make<TH1D>("fitcon", "fitcon", 1000, 0., 1.);
    momerr_ = dir.make<TH1D>("momerr", "momerr", 100, 0., 2.);
    t0err_ = dir.make<TH1D>("t0err", "t0err", 100, 0., 5.);
    t0_ = dir.make<TH1D>("t0", "t0", 339, 0., 1695.);
    d0_ = dir.make<TH1D>("d0", "d0", 200, -500., +500.);
    dout_ = dir.make<TH1D>("dout", "dout", 200, 0., +1000.);
    tanDip_ = dir.make<TH1D>("tanDip", "tanDip", 100, 0., 2.);
    momentum_ = dir.make<TH1D>("momentum", "momentum", 120, 100., 106.);
  }

  bool TrackCuts::accepted(const TrackSummary& trk, double weight) {
    hall_.fill(trk, weight);
    const auto& h = trk.states().at(0).helix();

    int numfailed=0, lastfailed = -1;

    if(trk.nactive() < cutnactive_) {
      ++numfailed;
      lastfailed = 0;
    }

    if(trk.fitcon() <= cutfitcon_) {
      ++numfailed;
      lastfailed = 1;
    }

    if(trk.states().at(0).momentumError() >= cutmomerr_) {
      ++numfailed;
      lastfailed = 2;
    }

    if(trk.t0Err() >= cutt0err_) {
      ++numfailed;
      lastfailed = 3;
    }

    if((h.d0() <= cutd0min_) || (cutd0max_ <= h.d0())) {
      ++numfailed;
      lastfailed = 4;
    }

    if((h.dOut() <= cutdoutmin_) || (cutdoutmax_ <= h.dOut())) {
      ++numfailed;
      lastfailed = 5;
    }

    const bool timeCutPassed =
      (cutt0min_ < cutt0max_) ?  // check for time wrapping
      ((cutt0min_ < trk.t0()) && (trk.t0() < cutt0max_)) : // no wrapping
      ((trk.t0() < cutt0max_) || (cutt0min_ < trk.t0()));  // time wrapped, discontinuous acceptance region

    if(!timeCutPassed) {
      ++numfailed;
      lastfailed = 6;
    }

    if((h.tanDip() <= cuttandipmin_) || (cuttandipmax_ <= h.tanDip())) {
      ++numfailed;
      lastfailed = 7;
    }

    const double pmag = trk.states().at(0).momentum().mag();
    if((pmag <= cutmommin_) || (cutmommax_ <= pmag)) {
      ++numfailed;
      lastfailed = 8;
    }

    if(numfailed == 0) {
      hfinal_.fill(trk, weight);

      nactive_->Fill(trk.nactive(), weight);
      fitcon_->Fill(trk.fitcon(), weight);
      momerr_->Fill(trk.states().at(0).momentumError(), weight);
      t0err_->Fill(trk.t0Err(), weight);
      d0_->Fill(h.d0(), weight);
      dout_->Fill(h.dOut(), weight);
      t0_->Fill(trk.t0(), weight);
      tanDip_->Fill(h.tanDip(), weight);
      momentum_->Fill(trk.states().at(0).momentum().mag(), weight);
    }
    else if(numfailed == 1) {
      switch(lastfailed) {
      case 0: nactive_->Fill(trk.nactive(), weight); break;
      case 1: fitcon_->Fill(trk.fitcon(), weight); break;
      case 2: momerr_->Fill(trk.states().at(0).momentumError(), weight); break;
      case 3: t0err_->Fill(trk.t0Err(), weight); break;
      case 4: d0_->Fill(h.d0(), weight); break;
      case 5: dout_->Fill(h.dOut(), weight); break;
      case 6: t0_->Fill(trk.t0(), weight); break;
      case 7: tanDip_->Fill(h.tanDip(), weight); break;
      case 8: momentum_->Fill(trk.states().at(0).momentum().mag(), weight); break;
      default: throw cet::exception("BUG")<<"TrackSummaryAnalyzerModule: unknown cut number";
      }
    }

    return (numfailed == 0);
  }

}
