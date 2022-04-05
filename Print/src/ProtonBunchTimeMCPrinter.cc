
#include "Offline/Print/inc/ProtonBunchTimeMCPrinter.hh"
#include "art/Framework/Principal/Provenance.h"
#include <string>
#include <iomanip>

void
mu2e::ProtonBunchTimeMCPrinter::Print(art::Event const& event,
                                std::ostream& os) {
  if(verbose()<1) return;
  if(tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector< art::Handle<ProtonBunchTimeMC> > vah = event.getMany<ProtonBunchTimeMC>();
    for (auto const & ah : vah) Print(ah);
  } else {
    // print requested instances
    for(const auto& tag : tags() ) {
      auto ih = event.getValidHandle<ProtonBunchTimeMC>(tag);
      Print(ih);
    }
  }
}

void
mu2e::ProtonBunchTimeMCPrinter::Print(const art::Handle<ProtonBunchTimeMC>& handle,
                                std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void
mu2e::ProtonBunchTimeMCPrinter::Print(const art::ValidHandle<ProtonBunchTimeMC>& handle,
                                std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}


void
mu2e::ProtonBunchTimeMCPrinter::Print(const mu2e::ProtonBunchTimeMC& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;

  os << std::setiosflags(std::ios::fixed | std::ios::right);

  os
    << " time: " << std::setw(8) << std::setprecision(2) << obj.pbtime_;
  os << std::endl;

}

void
mu2e::ProtonBunchTimeMCPrinter::PrintHeader(const std::string& tag, std::ostream& os) {
  if(verbose()<1) return;
  os << "\nProductPrint " << tag << "\n";
}


