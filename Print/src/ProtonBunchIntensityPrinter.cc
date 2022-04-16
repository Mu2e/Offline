
#include "Offline/Print/inc/ProtonBunchIntensityPrinter.hh"
#include "art/Framework/Principal/Provenance.h"
#include <iomanip>
#include <string>

void mu2e::ProtonBunchIntensityPrinter::Print(art::Event const& event,
                                              std::ostream& os) {
  if (verbose() < 1) return;
  if (tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector<art::Handle<ProtonBunchIntensity> > vah =
        event.getMany<ProtonBunchIntensity>();
    for (auto const& ah : vah) Print(ah);
  } else {
    // print requested instances
    for (const auto& tag : tags()) {
      auto ih = event.getValidHandle<ProtonBunchIntensity>(tag);
      Print(ih);
    }
  }
}

void mu2e::ProtonBunchIntensityPrinter::Print(
    const art::Handle<ProtonBunchIntensity>& handle, std::ostream& os) {
  if (verbose() < 1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back();  // remove trailing dot
  PrintHeader(tag, os);
  Print(*handle);
}

void mu2e::ProtonBunchIntensityPrinter::Print(
    const art::ValidHandle<ProtonBunchIntensity>& handle, std::ostream& os) {
  if (verbose() < 1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back();  // remove trailing dot
  PrintHeader(tag, os);
  Print(*handle);
}

void mu2e::ProtonBunchIntensityPrinter::Print(
    const mu2e::ProtonBunchIntensity& obj, int ind, std::ostream& os) {
  if (verbose() < 1) return;

  os << std::setiosflags(std::ios::fixed | std::ios::right);

  os << "  intensity: " << obj.intensity() << std::endl;
}

void mu2e::ProtonBunchIntensityPrinter::PrintHeader(const std::string& tag,
                                                    std::ostream& os) {
  if (verbose() < 1) return;
  os << "\nProductPrint " << tag << "\n";
}
