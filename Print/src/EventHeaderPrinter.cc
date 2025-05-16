
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"
#include "Offline/Mu2eUtilities/inc/EventHeaderFacade.hh"
#include "Offline/Print/inc/EventHeaderPrinter.hh"
#include "art/Framework/Principal/Provenance.h"
#include <iomanip>
#include <string>

void mu2e::EventHeaderPrinter::Print(art::Event const& event, std::ostream& os) {
  if (verbose() < 1) return;
  if (tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector<art::Handle<EventHeader> > vah = event.getMany<EventHeader>();
    for (auto const& ah : vah) Print(ah);
  } else {
    // print requested instances
    for (const auto& tag : tags()) {
      auto ih = event.getValidHandle<EventHeader>(tag);
      Print(ih);
    }
  }
}

void mu2e::EventHeaderPrinter::Print(const art::Handle<EventHeader>& handle,
                                  std::ostream& os) {
  if (verbose() < 1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back();  // remove trailing dot
  PrintHeader(tag, os);
  Print(*handle);
}

void mu2e::EventHeaderPrinter::Print(const art::ValidHandle<EventHeader>& handle,
                                  std::ostream& os) {
  if (verbose() < 1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back();  // remove trailing dot
  PrintHeader(tag, os);
  Print(*handle);
}

void mu2e::EventHeaderPrinter::Print(const mu2e::EventHeader& obj, int ind,
                                  std::ostream& os) {
  if (verbose() < 1) return;

  os << std::setiosflags(std::ios::fixed | std::ios::right);

  os << obj << std::endl;

  if (verbose() < 2) return;
  EventHeaderFacade ehf( obj, *(GlobalConstantsHandle<PhysicsParams>()) );
   os <<   "  Event Duration:       " << ehf.eventDuration() << " ns"
      << "\n  RF0 Offset Estimated: " << ehf.rf0OffsetEstimated() << " ns"
      << "\n  RF0 Offset Measured:  " << ehf.rf0OffsetMeasured() << " ns"
      << std::endl;
}

void mu2e::EventHeaderPrinter::PrintHeader(const std::string& tag,
                                        std::ostream& os) {
  if (verbose() < 1) return;
  os << "\nProductPrint " << tag << "\n";
}
