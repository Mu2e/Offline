#include "Offline/Print/inc/ProtonBunchIntensityPrinter.hh"
#include "art/Framework/Principal/Provenance.h"
#include <iomanip>
#include <string>

void mu2e::ProtonBunchIntensityPrinter::Print(art::Event const& event,
                                              std::ostream& os) {
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
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back();  // remove trailing dot
  PrintHeader(tag, os);
  Print(*handle);
}

void mu2e::ProtonBunchIntensityPrinter::Print(
    const art::ValidHandle<ProtonBunchIntensity>& handle, std::ostream& os) {
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back();  // remove trailing dot
  PrintHeader(tag, os);
  Print(*handle);
}

void mu2e::ProtonBunchIntensityPrinter::Print(
    const mu2e::ProtonBunchIntensity& obj, int ind, std::ostream& os) {
  if (verbose() >0 ){
    os << std::setiosflags(std::ios::fixed | std::ios::right);
    os << "  intensity: " << std::scientific << obj.intensity() << std::endl;
  }
  // summary
  nevts_++;
  nPOT_ += obj.intensity();
}

void mu2e::ProtonBunchIntensityPrinter::PrintHeader(const std::string& tag,
                                                    std::ostream& os) {
  if (verbose() < 1) return;
  os << "\nProductPrint " << tag << "\n";
}

void mu2e::ProtonBunchIntensityPrinter::PrintEndJob(std::ostream& os) {
  if(nevts_ > 0){
    os << "Processed " << nevts_ << " events for a total of " << std::scientific << std::setprecision(3) << nPOT_ << " POT" << std::endl;
  }
}
