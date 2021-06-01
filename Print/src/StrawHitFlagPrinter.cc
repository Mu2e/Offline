
#include "Print/inc/StrawHitFlagPrinter.hh"
#include "art/Framework/Principal/Provenance.h"
#include <string>
#include <iomanip>

void 
mu2e::StrawHitFlagPrinter::Print(art::Event const& event,
				std::ostream& os) {
  if(verbose()<1) return;
  if(tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector< art::Handle<StrawHitFlagCollection> > vah = event.getMany<StrawHitFlagCollection>();
    for (auto const & ah : vah) Print(ah);
  } else {
    // print requested instances
    for(const auto& tag : tags() ) {
      auto ih = event.getValidHandle<StrawHitFlagCollection>(tag);
      Print(ih);
    }
  }
}

void 
mu2e::StrawHitFlagPrinter::Print(const art::Handle<StrawHitFlagCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::StrawHitFlagPrinter::Print(const art::ValidHandle<StrawHitFlagCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::StrawHitFlagPrinter::Print(const StrawHitFlagCollection& coll, std::ostream& os) {
  if(verbose()<1) return;
  os << "StrawHitFlagCollection has " << coll.size() << " hits\n";
  if(verbose()==1) PrintListHeader();
  int i = 0;
  for(const auto& obj: coll) Print(obj, i++);
}

void 
mu2e::StrawHitFlagPrinter::Print(const mu2e::StrawHitFlag& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;

  os << std::setiosflags(std::ios::fixed | std::ios::right);
  if(ind>=0) os << std::setw(4) << ind;

  for(auto sn: obj.bitNames()) { 
    if(obj.hasAnyProperty(StrawHitFlag(sn.first))) 
      os << " " << sn.first;
  }
  os << std::endl;

}

void 
mu2e::StrawHitFlagPrinter::PrintHeader(const std::string& tag, std::ostream& os) {
  if(verbose()<1) return;
  os << "\nProductPrint " << tag << "\n";
}

void 
mu2e::StrawHitFlagPrinter::PrintListHeader(std::ostream& os) {
  if(verbose()<1) return;
  os << " ind    flags\n";

}

