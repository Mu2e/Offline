
#include "Print/inc/SimParticleTimeMapPrinter.hh"
#include "art/Framework/Principal/Provenance.h"
#include <string>
#include <iomanip>

void 
mu2e::SimParticleTimeMapPrinter::Print(art::Event const& event,
				std::ostream& os) {
  if(verbose()<1) return;
  if(tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector< art::Handle<SimParticleTimeMap> > vah = event.getMany<SimParticleTimeMap>();
    for (auto const & ah : vah) Print(ah);
  } else {
    // print requested instances
    for(const auto& tag : tags() ) {
      auto ih = event.getValidHandle<SimParticleTimeMap>(tag);
      Print(ih);
    }
  }
}

void 
mu2e::SimParticleTimeMapPrinter::Print(const art::Handle<SimParticleTimeMap>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::SimParticleTimeMapPrinter::Print(const art::ValidHandle<SimParticleTimeMap>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::SimParticleTimeMapPrinter::Print(const SimParticleTimeMap& coll, std::ostream& os) {
  if(verbose()<1) return;
  os << "SimParticleTimeMap has " << coll.size() << " times\n";
  if(verbose()!=1) return;

  os << std::setiosflags(std::ios::fixed | std::ios::right);

  PrintListHeader();

  int i = 0;
  long unsigned int pkey;
  for(const auto& obj: coll) {
    pkey = 0;
    art::Ptr<SimParticle> const& pptr = obj.first;
    if(pptr) pkey = pptr->id().asUint();
    os 
      << " " << std::setw(5) << i
      << " " << std::setw(7) << pkey
      << " " 
      << " " << std::setw(6) << std::setprecision(1) << obj.second
      << std::endl;
  }

}



void 
mu2e::SimParticleTimeMapPrinter::PrintHeader(const std::string& tag, std::ostream& os) {
  if(verbose()<1) return;
  os << "\nProductPrint " << tag << "\n";
}

void 
mu2e::SimParticleTimeMapPrinter::PrintListHeader(std::ostream& os) {
  if(verbose()<1) return;
  os << "  ind      id     time\n";

}

