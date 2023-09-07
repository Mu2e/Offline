#include "Offline/Print/inc/CosmicTrackSeedPrinter.hh"
#include "art/Framework/Principal/Provenance.h"
#include <iomanip>
#include <string>

void mu2e::CosmicTrackSeedPrinter::Print(art::Event const& event,
                                     std::ostream& os) {
  if (verbose() < 1) return;
  if (tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector<art::Handle<CosmicTrackSeedCollection> > vah =
        event.getMany<CosmicTrackSeedCollection>();
    for (auto const& ah : vah) Print(ah);
  } else {
    // print requested instances
    for (const auto& tag : tags()) {
      auto ih = event.getValidHandle<CosmicTrackSeedCollection>(tag);
      Print(ih);
    }
  }
}

void mu2e::CosmicTrackSeedPrinter::Print(
    const art::Handle<CosmicTrackSeedCollection>& handle, std::ostream& os) {
  if (verbose() < 1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back();  // remove trailing dot
  PrintHeader(tag, os);
  Print(*handle);
}

void mu2e::CosmicTrackSeedPrinter::Print(
    const art::ValidHandle<CosmicTrackSeedCollection>& handle, std::ostream& os) {
  if (verbose() < 1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back();  // remove trailing dot
  PrintHeader(tag, os);
  Print(*handle);
}

void mu2e::CosmicTrackSeedPrinter::Print(const CosmicTrackSeedCollection& coll,
                                     std::ostream& os) {
  if (verbose() < 1) return;
  os << "CosmicTrackSeedCollection has " << coll.size() << " CosmicTrack Seeds\n";
  if (verbose() >= 1) PrintListHeader();
  int i = 0;
  for (const auto& obj : coll) Print(obj, i++);
}

void mu2e::CosmicTrackSeedPrinter::Print(const art::Ptr<CosmicTrackSeed>& obj, int ind,
                                     std::ostream& os) {
  if (verbose() < 1) return;
  Print(*obj, ind);
}

void mu2e::CosmicTrackSeedPrinter::Print(const mu2e::CosmicTrackSeed& obj, int ind,
                                     std::ostream& os) {
  if (verbose() < 1) return;

  os << std::setiosflags(std::ios::fixed | std::ios::right);
  if (ind >= 0) os << std::setw(4) << ind;

  os << " " << std::setw(5) << obj.hits().size()
     << " " << std::setw(5) << obj.status()
     << " " << std::setw(8) << std::setprecision(3) << obj.track().MinuitParams.A0
     << " " << std::setw(8) << std::setprecision(3) << obj.track().MinuitParams.A1
     << " " << std::setw(8) << std::setprecision(3) << obj.track().MinuitParams.B0
     << " " << std::setw(8) << std::setprecision(3) << obj.track().MinuitParams.B1
     << " " << std::setw(8) << std::setprecision(3) << obj.t0().t0() << std::endl;

  if(verbose() > 2){
    for(auto const& hit : obj.hits()) {
      os << hit << std::endl;
    }
  }
}

void mu2e::CosmicTrackSeedPrinter::PrintHeader(const std::string& tag,
                                           std::ostream& os) {
  if (verbose() < 1) return;
  os << "\nProductPrint " << tag << "\n";
}

void mu2e::CosmicTrackSeedPrinter::PrintListHeader(std::ostream& os) {
  if (verbose() < 1) return;
  os << " ind   nhits   status  A0      A1      B0      B1      t0 \n";
}
