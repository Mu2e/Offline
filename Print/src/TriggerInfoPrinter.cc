
#include "Offline/Print/inc/TriggerInfoPrinter.hh"
#include "art/Framework/Principal/Provenance.h"
#include <iomanip>
#include <string>
#include <vector>

void mu2e::TriggerInfoPrinter::Print(art::Event const& event,
                                     std::ostream& os) {
  if (verbose() < 1) return;
  if (tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector<art::Handle<mu2e::TriggerInfo> > vah =
        event.getMany<mu2e::TriggerInfo>();
    for (auto const& ah : vah) Print(ah);
  } else {
    // print requested instances
    for (const auto& tag : tags()) {
      auto ih = event.getValidHandle<mu2e::TriggerInfo>(tag);
      Print(ih);
    }
  }
}

void mu2e::TriggerInfoPrinter::Print(
    const art::Handle<mu2e::TriggerInfo>& handle, std::ostream& os) {
  if (verbose() < 1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back();  // remove trailing dot
  PrintHeader(tag, os);
  Print(*handle);
}

void mu2e::TriggerInfoPrinter::Print(
    const art::ValidHandle<mu2e::TriggerInfo>& handle, std::ostream& os) {
  if (verbose() < 1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back();  // remove trailing dot
  PrintHeader(tag, os);
  Print(*handle);
}

void mu2e::TriggerInfoPrinter::Print(const mu2e::TriggerInfo& obj, int ind,
                                     std::ostream& os) {
  if (verbose() < 1) return;

  os << std::setiosflags(std::ios::fixed | std::ios::right);
  std::string avail;
  avail = "";
  if (obj.caloClusters().size() == 0) avail = ", some unavailable";
  os << std::setw(4) << obj.caloClusters().size() << " caloClusters" << avail
     << std::endl;
  avail = "";
  if (obj.tracks().size() == 0) avail = ", some unavailable";
  os << std::setw(4) << obj.tracks().size() << " tracks" << avail << std::endl;
  avail = "";
  if (obj.helixes().size() ==0) avail = ", some unavailable";
  os << std::setw(4) << obj.helixes().size() << " helixes" << avail
     << std::endl;
  avail = "";
  if (obj.hitClusters().size() == 0) avail = ", some unavailable";
  os << std::setw(4) << obj.hitClusters().size() << " hitClusters" << avail
     << std::endl;
  avail = "";
  if (obj.caloTrigSeeds().size() == 0) avail = ", some unavailable";
  os << std::setw(4) << obj.caloTrigSeeds().size() << " caloTrigSeeds" << avail
     << std::endl;
  avail = "";
  if (obj.cosmics().size() == 0) avail = ", some unavailable";
  os << std::setw(4) << obj.cosmics().size() << " cosmics" << avail
     << std::endl;
}

void mu2e::TriggerInfoPrinter::PrintHeader(const std::string& tag,
                                           std::ostream& os) {
  if (verbose() < 1) return;
  os << "\nProductPrint " << tag << "\n";
}
