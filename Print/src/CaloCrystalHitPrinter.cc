
#include "Print/inc/CaloCrystalHitPrinter.hh"
#include "art/Framework/Principal/Provenance.h"
#include <string>
#include <iomanip>

void 
mu2e::CaloCrystalHitPrinter::Print(art::Event const& event,
				std::ostream& os) {
  if(verbose()<1) return;
  if(tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector< art::Handle<CaloCrystalHitCollection> > vah;
    event.getManyByType(vah);
    for (auto const & ah : vah) Print(ah);
  } else {
    // print requested instances
    for(const auto& tag : tags() ) {
      auto ih = event.getValidHandle<CaloCrystalHitCollection>(tag);
      Print(ih);
    }
  }
}

void 
mu2e::CaloCrystalHitPrinter::Print(const art::Handle<CaloCrystalHitCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::CaloCrystalHitPrinter::Print(const art::ValidHandle<CaloCrystalHitCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::CaloCrystalHitPrinter::Print(const CaloCrystalHitCollection& coll, std::ostream& os) {
  if(verbose()<1) return;
  os << "CaloCrystalHitCollection has " << coll.size() << " hits\n";
  if(verbose()==1) PrintListHeader();
  int i = 0;
  for(const auto& obj: coll) Print(obj, i++);
}

void 
mu2e::CaloCrystalHitPrinter::Print(const art::Ptr<CaloCrystalHit>& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;
  Print(*obj,ind);
}

void 
mu2e::CaloCrystalHitPrinter::Print(const mu2e::CaloCrystalHit& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;

  if( obj.energyDep() < _eCut ) return;

  os << std::setiosflags(std::ios::fixed | std::ios::right);
  if(ind>=0) os << std::setw(4) << ind;

  os 
    << " " << std::setw(5) << obj.id()
    << " " 
    << " " << std::setw(8) << std::setprecision(1) << obj.time()
    << " " << std::setw(8) << std::setprecision(1) << obj.energyDep()
    << std::endl;

}

void 
mu2e::CaloCrystalHitPrinter::PrintHeader(const std::string& tag, std::ostream& os) {
  if(verbose()<1) return;
  os << "\nProductPrint " << tag << "\n";
}

void 
mu2e::CaloCrystalHitPrinter::PrintListHeader(std::ostream& os) {
  if(verbose()<1) return;
  os << "ind     id     time     energy    eTot   nROIds\n";

}


