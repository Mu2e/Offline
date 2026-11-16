#include "Offline/Print/inc/RawEventHeaderPrinter.hh"
#include "art/Framework/Principal/Provenance.h"
#include <iomanip>
#include <string>

void mu2e::RawEventHeaderPrinter::Print(art::Event const& event,
                                        std::ostream& os) {
  if (verbose() < 1) return;
  if (tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector<art::Handle<artdaq::detail::RawEventHeader> > vah =
        event.getMany<artdaq::detail::RawEventHeader>();
    for (auto const& ah : vah) Print(ah);
  } else {
    // print requested instances
    for (const auto& tag : tags()) {
      auto ih = event.getValidHandle<artdaq::detail::RawEventHeader>(tag);
      Print(ih);
    }
  }
}

void mu2e::RawEventHeaderPrinter::Print(
    const art::Handle<artdaq::detail::RawEventHeader>& handle,
    std::ostream& os) {
  if (verbose() < 1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back();  // remove trailing dot
  PrintHeader(tag, os);
  Print(*handle, os);
}

void mu2e::RawEventHeaderPrinter::Print(
    const art::ValidHandle<artdaq::detail::RawEventHeader>& handle,
    std::ostream& os) {
  if (verbose() < 1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back();  // remove trailing dot
  PrintHeader(tag, os);
  Print(*handle, os);
}

void mu2e::RawEventHeaderPrinter::Print(
    const artdaq::detail::RawEventHeader& obj, std::ostream& os) {
  if (verbose() < 1) return;

  os << obj << std::endl;
}

void mu2e::RawEventHeaderPrinter::PrintHeader(const std::string& tag,
                                              std::ostream& os) {
  if (verbose() < 1) return;
  os << "\nProductPrint " << tag << "\n";
}
