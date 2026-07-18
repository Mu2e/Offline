#include "Offline/Print/inc/EventHeadersPrinter.hh"
#include "art/Framework/Principal/Provenance.h"
#include <iomanip>
#include <string>

void mu2e::EventHeadersPrinter::Print(art::Event const& event,
                                         std::ostream& os) {
  if (verbose() < 1) return;
  if (tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector<art::Handle<EventHeaders> > vah =
        event.getMany<EventHeaders>();
    for (auto const& ah : vah) Print(ah);
  } else {
    // print requested instances
    for (const auto& tag : tags()) {
      auto ih = event.getValidHandle<EventHeaders>(tag);
      Print(ih);
    }
  }
}

void mu2e::EventHeadersPrinter::Print(
    const art::Handle<EventHeaders>& handle, std::ostream& os) {
  if (verbose() < 1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back();  // remove trailing dot
  PrintHeader(tag, os);
  Print(*handle);
}

void mu2e::EventHeadersPrinter::Print(
    const art::ValidHandle<EventHeaders>& handle, std::ostream& os) {
  if (verbose() < 1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back();  // remove trailing dot
  PrintHeader(tag, os);
  Print(*handle);
}

void mu2e::EventHeadersPrinter::Print(const mu2e::EventHeaders& obj,
                                         int ind, std::ostream& os) {
  if (verbose() < 1) return;
  for (const auto& sobj : obj ) {
    os << sobj << std::endl;
  }

}

void mu2e::EventHeadersPrinter::PrintHeader(const std::string& tag,
                                               std::ostream& os) {
  if (verbose() < 1) return;
  os << "\nProductPrint " << tag << "\n";
}
