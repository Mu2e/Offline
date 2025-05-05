
#include "Offline/Print/inc/KalSeedPrinter.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "art/Framework/Principal/Provenance.h"
#include <string>

void mu2e::KalSeedPrinter::Print(art::Event const& event, std::ostream& os) {
  if (verbose() < 1) return;
  if (tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector<art::Handle<KalSeedCollection> > vah =
      event.getMany<KalSeedCollection>();
    for (auto const& ah : vah) Print(ah);
  } else {
    // print requested instances
    for (const auto& tag : tags()) {
      auto ih = event.getValidHandle<KalSeedCollection>(tag);
      Print(ih);
    }
  }
}

void mu2e::KalSeedPrinter::Print(const art::Handle<KalSeedCollection>& handle,
    std::ostream& os) {
  if (verbose() < 1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back();  // remove trailing dot
  PrintHeader(tag, os);
  Print(*handle);
}

void mu2e::KalSeedPrinter::Print(
    const art::ValidHandle<KalSeedCollection>& handle, std::ostream& os) {
  if (verbose() < 1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back();  // remove trailing dot
  PrintHeader(tag, os);
  Print(*handle);
}

void mu2e::KalSeedPrinter::Print(const KalSeedCollection& coll,
    std::ostream& os) {
  if (verbose() < 1) return;
  os << "KalSeedCollection has " << coll.size() << " tracks\n";
  if (verbose() == 1) PrintListHeader();
  int i = 0;
  for (const auto& obj : coll) Print(obj, i++);
}

void mu2e::KalSeedPrinter::Print(const art::Ptr<KalSeed>& obj, int ind,
    std::ostream& os) {
  if (verbose() < 1) return;
  Print(*obj, ind);
}

void mu2e::KalSeedPrinter::Print(const mu2e::KalSeed& obj, int ind,
    std::ostream& os) {
  if (verbose() < 1) return;

  os << std::setiosflags(std::ios::fixed | std::ios::right);
  if (ind >= 0 && verbose() == 1) os << std::setw(4) << ind;

  KalSegment seg;  // this will be filled with 0's and -1's
  // use the first segment, at the front of the tracker
  if (obj.segments().size() > 0) seg = obj.segments()[0];

  const mu2e::HelixVal& hh = seg.helix();

  if (verbose() == 1) {
    os << " " << std::setw(5) << obj.status().hex() << " " << std::setw(8)
      << std::setprecision(3) << obj.fitConsistency() << " " << std::setw(8)
      << std::setprecision(3) << seg.mom() << " " << std::setw(6)
      << std::setprecision(3) << seg.momerr() << " " << std::setw(8)
      << std::setprecision(4) << hh.tanDip() << " " << std::setw(7)
      << std::setprecision(1) << hh.d0() << " " << std::setw(7)
      << std::setprecision(5) << hh.omega() << " " << std::setw(7)
      << std::setprecision(1) << obj.t0().t0() << " " << std::setw(7)
      << std::setprecision(4) << obj.hits().size() << std::endl;
  } else if (verbose() >= 2) {
    auto const& ptable = GlobalConstantsHandle<ParticleDataList>();

    os << " fitStatus: " << std::setw(3) << obj.status() << "\n";
    os << " part: " << ptable->particle(obj.particle()).name()
       << "  t0: " << std::setw(7)
       << std::setprecision(1) << obj.t0().t0()
       << " chi2: " << std::setw(7) << std::setprecision(2) << obj.chisquared()
       << "  fitcon: " << std::setw(7) << std::setprecision(3)
       << obj.fitConsistency() << "  nhits: " << std::setw(3)
       << obj.hits().size()
       << "  calo: " << (obj.caloCluster().isNull() ? "no" : "yes") << "\n";
    os << " intersections: \n";
    for (auto const& inter : obj.intersections()) {
      os << " sid " << inter.surfaceId() << " time " << inter.time() << " P " << inter.momentum3() << " dP " << inter.dMom() << "\n";
    }
    if (verbose() >= 3) {
      os << " segments: \n";
      for (auto const& ss : obj.segments()) {
        const mu2e::HelixVal& h = ss.helix();
        const mu2e::HelixCov& w = ss.covar();
        CLHEP::HepSymMatrix c;
        w.symMatrix(c);  // fills c withthe cov

        os << " p: " << std::setw(8) << std::setprecision(3) << ss.mom()
          << "  +/- " << std::setw(6) << std::setprecision(3) << ss.momerr()
          << "    fmin: " << std::setw(8) << std::setprecision(6) << ss.fmin()
          << "  fmax: " << std::setw(6) << std::setprecision(1) << ss.fmax()
          << "\n";
        if (verbose() == 3) {
          os << "   d0: " << std::setw(5) << std::setprecision(1) << h.d0()
            << "  phi0: " << std::setw(6) << std::setprecision(3) << h.phi0()
            << "  omega: " << std::setw(8) << std::setprecision(6) << h.omega()
            << "  z0: " << std::setw(6) << std::setprecision(1) << h.z0()
            << "  tanDip: " << std::setw(7) << std::setprecision(3) << h.tanDip()
            << "\n";
        } else {
          os << "     d0: " << std::setw(8) << std::setprecision(1) << h.d0()
            << " +/- " << std::setw(9) << std::setprecision(2) << sqrt(c(1, 1))
            << "\n";
          os << "   phi0: " << std::setw(8) << std::setprecision(4) << hh.phi0()
            << " +/- " << std::setw(9) << std::setprecision(5) << sqrt(c(2, 2))
            << "\n";
          os << "  omega: " << std::setw(8) << std::setprecision(6) << hh.omega()
            << " +/- " << std::setw(9) << std::setprecision(7) << sqrt(c(3, 3))
            << "\n";
          os << "     z0: " << std::setw(8) << std::setprecision(1) << hh.z0()
            << " +/- " << std::setw(9) << std::setprecision(2) << sqrt(c(4, 4))
            << "\n";
          os << " tanDip: " << std::setw(8) << std::setprecision(4) << hh.tanDip()
            << " +/- " << std::setw(9) << std::setprecision(5) << sqrt(c(5, 5))
            << "\n";
        }
        if (verbose() >= 4) {
          os << "  Helix covariance:\n";
          PrintMatrix(c, os, 0);

          os << "  Helix correlations:\n";
          PrintMatrix(c, os, 1);
        }
      }
    }
  }
}

void mu2e::KalSeedPrinter::PrintHeader(const std::string& tag,
    std::ostream& os) {
  if (verbose() < 1) return;
  os << "\nProductPrint " << tag << "\n";
}

void mu2e::KalSeedPrinter::PrintListHeader(std::ostream& os) {
  if (verbose() < 1) return;
  os << "ind  status   fitcon    p      pErr   tanDip     d0     omega    t0   "
    " nhits\n";
}
