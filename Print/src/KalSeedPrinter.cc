
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
  double t0;
  auto seg = obj.t0Segment(t0);

  if (verbose() >= 1) {
    std::string fittype("Unknown");
    if(obj.loopHelixFit())fittype = "LoopHelix";
    if(obj.centralHelixFit())fittype = "CentralHelix";
    if(obj.kinematicLineFit())fittype = "KinematicLine";
    unsigned nactive(0);
    for(auto const& hit : obj.hits())if(hit._ambig>WireHitState::inactive)nactive++;
    os << " " << std::setprecision(2) << obj.particle() << " " << std::setw(2)
      << std::setprecision(2) << fittype << " " << std::setw(2)
      << std::setprecision(3) << obj.fitConsistency() << " " << std::setw(5)
      << std::setprecision(1) << (obj.caloCluster().isNull() ? "no" : "yes") << std::setw(3)
      << std::setprecision(4) << nactive << " " << std::setw(3)
      << std::setprecision(4) << obj.straws().size() << " " << std::setw(3)
      << std::setprecision(4) << obj.segments().size() << " " << std::setw(3)
      << std::setprecision(4) << obj.intersections().size() << " " << std::setw(3);
    if(seg != obj.segments().end()){
      os << std::setprecision(3) << seg->mom() << " " << std::setw(6)
      << std::setprecision(3) << seg->momerr() << " " << std::setw(8)
      << std::setprecision(3) << cos(seg->momentum3().Theta()) << " " << std::setw(7)
      << std::setprecision(1) << t0 << " " << std::setw(7);
    }
    os << "\n";

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
  }
}

void mu2e::KalSeedPrinter::PrintHeader(const std::string& tag,
    std::ostream& os) {
  if (verbose() < 1) return;
  os << "\nProductPrint " << tag << "\n";
}

void mu2e::KalSeedPrinter::PrintListHeader(std::ostream& os) {
  if (verbose() < 1) return;
  os << "ind pdg fittype  fitcon calo? nhit nst nseg nint p   pErr cos(theta)  t0 \n";
}
