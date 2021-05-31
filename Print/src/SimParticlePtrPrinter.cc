
#include "Print/inc/SimParticlePtrPrinter.hh"
#include "art/Framework/Principal/Provenance.h"
#include <string>
#include <iomanip>

void 
mu2e::SimParticlePtrPrinter::Print(art::Event const& event,
				std::ostream& os) {
  if(verbose()<1) return;
  if(tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector< art::Handle<SimParticlePtrCollection> > vah = event.getMany<SimParticlePtrCollection>();
    for (auto const & ah : vah) Print(ah);
  } else {
    // print requested instances
    for(const auto& tag : tags() ) {
      auto ih = event.getValidHandle<SimParticlePtrCollection>(tag);
      Print(ih);
    }
  }
}

void 
mu2e::SimParticlePtrPrinter::Print(const art::Handle<SimParticlePtrCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::SimParticlePtrPrinter::Print(const art::ValidHandle<SimParticlePtrCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::SimParticlePtrPrinter::Print(const SimParticlePtrCollection& coll, std::ostream& os) {
  if(verbose()<1) return;
  os << "SimParticlePtrCollection has " << coll.size() << " pointers\n";
  if(verbose()==1) PrintListHeader();
  int i = 0;
  for(const auto& obj: coll) Print(obj, i++);
}

void 
mu2e::SimParticlePtrPrinter::Print(const art::Ptr<SimParticle>& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;

  os << std::setiosflags(std::ios::fixed | std::ios::right);
  if(ind>=0) os << std::setw(4) << ind;

  if(verbose()==1) {
    os 
      << " " << std::setw(12) << obj.id()
      << " " << std::setw(12) << obj.key()
      << std::endl;

  } // end if verbose

}

void 
mu2e::SimParticlePtrPrinter::PrintHeader(const std::string& tag, std::ostream& os) {
  if(verbose()<1) return;
  os << "\nProductPrint " << tag << "\n";
}

void 
mu2e::SimParticlePtrPrinter::PrintListHeader(std::ostream& os) {
  if(verbose()<1) return;
  os << "ind     ProductID        key\n";
}

