#include "Offline/Print/inc/HelixSeedPrinter.hh"
#include "art/Framework/Principal/Provenance.h"
#include <iomanip>
#include <string>

void mu2e::HelixSeedPrinter::Print(art::Event const& event,
                                     std::ostream& os) {
  if (verbose() < 1) return;
  if (tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector<art::Handle<HelixSeedCollection> > vah =
        event.getMany<HelixSeedCollection>();
    for (auto const& ah : vah) Print(ah);
  } else {
    // print requested instances
    for (const auto& tag : tags()) {
      auto ih = event.getValidHandle<HelixSeedCollection>(tag);
      Print(ih);
    }
  }
}

void mu2e::HelixSeedPrinter::Print(
    const art::Handle<HelixSeedCollection>& handle, std::ostream& os) {
  if (verbose() < 1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back();  // remove trailing dot
  PrintHeader(tag, os);
  Print(*handle);
}

void mu2e::HelixSeedPrinter::Print(
    const art::ValidHandle<HelixSeedCollection>& handle, std::ostream& os) {
  if (verbose() < 1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back();  // remove trailing dot
  PrintHeader(tag, os);
  Print(*handle);
}

void mu2e::HelixSeedPrinter::Print(const HelixSeedCollection& coll,
                                     std::ostream& os) {
  if (verbose() < 1) return;
  os << "HelixSeedCollection has " << coll.size() << " Helix Seeds\n";
  if (verbose() >= 1) PrintListHeader();
  int i = 0;
  for (const auto& obj : coll) Print(obj, i++);
}

void mu2e::HelixSeedPrinter::Print(const art::Ptr<HelixSeed>& obj, int ind,
                                     std::ostream& os) {
  if (verbose() < 1) return;
  Print(*obj, ind);
}

void mu2e::HelixSeedPrinter::Print(const mu2e::HelixSeed& obj, int ind,
                                     std::ostream& os) {
  if (verbose() < 1) return;

  os << std::setiosflags(std::ios::fixed | std::ios::right);
  if (ind >= 0) os << std::setw(4) << ind;

  os << " " << std::setw(5) << obj.hits().size()
     << " " << std::setw(5) << obj.status()
     << " " << std::setw(8) << std::setprecision(3) << obj.helix().centerx()
     << " " << std::setw(8) << std::setprecision(3) << obj.helix().centery()
     << " " << std::setw(8) << std::setprecision(3) << obj.helix().radius()
     << " " << std::setw(8) << std::setprecision(3) << obj.helix().lambda()
     << " " << std::setw(8) << std::setprecision(3) << obj.helix().fz0()
     << std::setw(7) << std::setprecision(1) << obj.t0().t0() << std::endl;
}

void mu2e::HelixSeedPrinter::PrintHeader(const std::string& tag,
                                           std::ostream& os) {
  if (verbose() < 1) return;
  os << "\nProductPrint " << tag << "\n";
}

void mu2e::HelixSeedPrinter::PrintListHeader(std::ostream& os) {
  if (verbose() < 1) return;
  os << " ind   nhits   status  centerx  centery  radius  lambda  fz0   t0 \n";
}
