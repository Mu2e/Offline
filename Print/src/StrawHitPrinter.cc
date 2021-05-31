
#include "Print/inc/StrawHitPrinter.hh"
#include "art/Framework/Principal/Provenance.h"
#include <string>
#include <iomanip>

void 
mu2e::StrawHitPrinter::Print(art::Event const& event,
				std::ostream& os) {
  if(verbose()<1) return;
  if(tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector< art::Handle<StrawHitCollection> > vah = event.getMany<StrawHitCollection>();
    for (auto const & ah : vah) Print(ah);
  } else {
    // print requested instances
    for(const auto& tag : tags() ) {
      auto ih = event.getValidHandle<StrawHitCollection>(tag);
      Print(ih);
    }
  }
}

void 
mu2e::StrawHitPrinter::Print(const art::Handle<StrawHitCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::StrawHitPrinter::Print(const art::ValidHandle<StrawHitCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::StrawHitPrinter::Print(const StrawHitCollection& coll, std::ostream& os) {
  if(verbose()<1) return;
  os << "StrawHitCollection has " << coll.size() << " hits\n";
  if(verbose()==1) PrintListHeader();
  int i = 0;
  for(const auto& obj: coll) Print(obj, i++);
}

void 
mu2e::StrawHitPrinter::Print(const art::Ptr<StrawHit>& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;
  Print(*obj,ind);
}

void 
mu2e::StrawHitPrinter::Print(const mu2e::StrawHit& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;

  if( obj.energyDep() < _eCut ) return;

  os << std::setiosflags(std::ios::fixed | std::ios::right);
  if(ind>=0) os << std::setw(4) << ind;

  os 
    << " " << std::setw(5) << obj.strawId().asUint16()
    << " " 
    << " " << std::setw(8) << std::setprecision(1) << obj.time()
    << " " << std::setw(8) << std::setprecision(3) << obj.dt()
    << " " << std::setw(10) << std::setprecision(6) << obj.energyDep()
    << std::endl;

}

void 
mu2e::StrawHitPrinter::PrintHeader(const std::string& tag, std::ostream& os) {
  if(verbose()<1) return;
  os << "\nProductPrint " << tag << "\n";
}

void 
mu2e::StrawHitPrinter::PrintListHeader(std::ostream& os) {
  if(verbose()<1) return;
  os << " ind StrwInd   time       dt      eDep\n";

}

