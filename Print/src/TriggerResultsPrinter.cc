
#include "Offline/Print/inc/TriggerResultsPrinter.hh"
#include "Offline/Mu2eUtilities/inc/TriggerResultsNavigator.hh"
#include "art/Framework/Principal/Provenance.h"
#include <iomanip>
#include <string>
#include <vector>

void mu2e::TriggerResultsPrinter::Print(art::Event const& event,
                                        std::ostream& os) {
  if (verbose() < 1) return;
  if (tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector<art::Handle<art::TriggerResults> > vah =
        event.getMany<art::TriggerResults>();
    for (auto const& ah : vah) Print(ah);
  } else {
    // print requested instances
    for (const auto& tag : tags()) {
      auto ih = event.getValidHandle<art::TriggerResults>(tag);
      Print(ih);
    }
  }
}

void mu2e::TriggerResultsPrinter::Print(
    const art::Handle<art::TriggerResults>& handle, std::ostream& os) {
  if (verbose() < 1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back();  // remove trailing dot
  PrintHeader(tag, os);
  Print(*handle);
}

void mu2e::TriggerResultsPrinter::Print(
    const art::ValidHandle<art::TriggerResults>& handle, std::ostream& os) {
  if (verbose() < 1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back();  // remove trailing dot
  PrintHeader(tag, os);
  Print(*handle);
}

void mu2e::TriggerResultsPrinter::Print(const art::TriggerResults& obj, int ind,
                                        std::ostream& os) {
  if (verbose() < 1) return;

  os << std::setiosflags(std::ios::fixed | std::ios::right);
  std::vector<std::string> text{"n", "1", "0", "e"};

  TriggerResultsNavigator trigNavig(&obj);
  for (unsigned int i = 0; i < trigNavig.getTrigPaths().size(); ++i) {
    std::string path = trigNavig.getTrigPathNameByIndex(i);
    size_t pathID = trigNavig.getTrigBitByName(path);
    os << "  " << (trigNavig.accepted(path) ? "pass" : "fail") << "  " << pathID
       << "  " << path << std::endl;
  }
}

void mu2e::TriggerResultsPrinter::PrintHeader(const std::string& tag,
                                              std::ostream& os) {
  if (verbose() < 1) return;
  os << "\nProductPrint " << tag << "\n";
}
