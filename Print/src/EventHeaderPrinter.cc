#include "Offline/Print/inc/EventHeaderPrinter.hh"
#include "art/Framework/Principal/Provenance.h"
#include <iomanip>
#include <string>

void mu2e::EventHeaderPrinter::Print(art::Event const& event,
                                         std::ostream& os) {
  if (verbose() < 1) return;
  if (tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector<art::Handle<EventHeader> > vah =
        event.getMany<EventHeader>();
    for (auto const& ah : vah) Print(ah, os);
  } else {
    // print requested instances
    for (const auto& tag : tags()) {
      auto ih = event.getValidHandle<EventHeader>(tag);
      Print(ih, os);
    }
  }
}

void mu2e::EventHeaderPrinter::Print(
    const art::Handle<EventHeader>& handle, std::ostream& os) {
  if (verbose() < 1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back();  // remove trailing dot
  PrintHeader(tag, os);
  Print(*handle, os);
}

void mu2e::EventHeaderPrinter::Print(
    const art::ValidHandle<EventHeader>& handle, std::ostream& os) {
  if (verbose() < 1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back();  // remove trailing dot
  PrintHeader(tag, os);
  Print(*handle, os);
}

void mu2e::EventHeaderPrinter::Print(const mu2e::EventHeader& obj,
                                         std::ostream& os) {
  if (verbose() < 1) return;

  os << obj << std::endl;
}

void mu2e::EventHeaderPrinter::PrintHeader(const std::string& tag,
                                               std::ostream& os) {
  if (verbose() < 1) return;
  os << "\nProductPrint " << tag << "\n";
}
