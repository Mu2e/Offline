#include "Print/inc/PrimaryParticlePrinter.hh"
#include "Print/inc/GenParticlePrinter.hh"
#include "Print/inc/SimParticlePrinter.hh"
#include "art/Framework/Principal/Provenance.h"
#include <string>
#include <iomanip>

void 
mu2e::PrimaryParticlePrinter::Print(art::Event const& event,
				std::ostream& os) {
  if(verbose()<1) return;
  if(tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector< art::Handle<PrimaryParticle> > vah = event.getMany<PrimaryParticle>();
    for (auto const & ah : vah) Print(ah);
  } else {
    // print requested instances
    for(const auto& tag : tags() ) {
      auto ih = event.getValidHandle<PrimaryParticle>(tag);
      Print(ih);
    }
  }
}

void 
mu2e::PrimaryParticlePrinter::Print(const art::Handle<PrimaryParticle>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::PrimaryParticlePrinter::Print(const art::ValidHandle<PrimaryParticle>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}


void 
mu2e::PrimaryParticlePrinter::Print(const mu2e::PrimaryParticle& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;
  os << std::setiosflags(std::ios::fixed | std::ios::right);
  os << "PrimaryParticle GenParticle:" << std::endl;
  gprint_.PrintListHeader(os);
  gprint_.Print(obj.primary(),-1,os);
  os << "PrimaryParticle has "<< obj.primarySimParticles().size() << " SimParticles:" << std::endl;
  sprint_.PrintListHeader(os);
  for(size_t isp = 0; isp < obj.primarySimParticles().size(); isp++){
    auto const& spp = obj.primarySimParticles()[isp];
    if(spp.isAvailable()) sprint_.Print(spp,isp,os);
  }
}

void 
mu2e::PrimaryParticlePrinter::PrintHeader(const std::string& tag, std::ostream& os) {
  if(verbose()<1) return;
  os << "\nProductPrint " << tag << "\n";
}

