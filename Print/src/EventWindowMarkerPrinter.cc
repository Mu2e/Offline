
#include "Offline/Print/inc/EventWindowMarkerPrinter.hh"
#include "art/Framework/Principal/Provenance.h"
#include <iomanip>
#include <string>

void mu2e::EventWindowMarkerPrinter::Print(art::Event const& event,
                                           std::ostream& os) {
  if (verbose() < 1) return;
  if (tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector<art::Handle<EventWindowMarker> > vah =
        event.getMany<EventWindowMarker>();
    for (auto const& ah : vah) Print(ah);
  } else {
    // print requested instances
    for (const auto& tag : tags()) {
      auto ih = event.getValidHandle<EventWindowMarker>(tag);
      Print(ih);
    }
  }
}

void mu2e::EventWindowMarkerPrinter::Print(
    const art::Handle<EventWindowMarker>& handle, std::ostream& os) {
  if (verbose() < 1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back();  // remove trailing dot
  PrintHeader(tag, os);
  Print(*handle);
}

void mu2e::EventWindowMarkerPrinter::Print(
    const art::ValidHandle<EventWindowMarker>& handle, std::ostream& os) {
  if (verbose() < 1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back();  // remove trailing dot
  PrintHeader(tag, os);
  Print(*handle);
}

void mu2e::EventWindowMarkerPrinter::Print(const mu2e::EventWindowMarker& obj,
                                           int ind, std::ostream& os) {
  if (verbose() < 1) return;

  os << std::setiosflags(std::ios::fixed | std::ios::right);

  os << " spillType: " << std::setw(3) << obj.spillType()
     << "   length: " << std::setw(8) << std::setprecision(2)
     << obj.eventLength();
  os << std::endl;
}

void mu2e::EventWindowMarkerPrinter::PrintHeader(const std::string& tag,
                                                 std::ostream& os) {
  if (verbose() < 1) return;
  os << "\nProductPrint " << tag << "\n";
}
