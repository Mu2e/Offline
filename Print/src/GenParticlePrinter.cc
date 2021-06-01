
#include "Print/inc/GenParticlePrinter.hh"
#include "art/Framework/Principal/Provenance.h"
#include <string>
#include <iomanip>

void 
mu2e::GenParticlePrinter::Print(art::Event const& event,
				std::ostream& os) {
  if(verbose()<1) return;
  if(tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector< art::Handle<GenParticleCollection> > vah = event.getMany<GenParticleCollection>();
    for (auto const & ah : vah) Print(ah);
  } else {
    // print requested instances
    for(const auto& tag : tags() ) {
      auto ih = event.getValidHandle<GenParticleCollection>(tag);
      Print(ih);
    }
  }
}

void 
mu2e::GenParticlePrinter::Print(const art::Handle<GenParticleCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::GenParticlePrinter::Print(const art::ValidHandle<GenParticleCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::GenParticlePrinter::Print(const GenParticleCollection& coll, std::ostream& os) {
  if(verbose()<1) return;
  os << "GenParticleCollection has " << coll.size() << " particles\n";
  if(verbose()==1) PrintListHeader();
  int i = 0;
  for(const auto& obj: coll) Print(obj, i++);
}

void 
mu2e::GenParticlePrinter::Print(const art::Ptr<GenParticle>& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;
  Print(*obj,ind);
}

void 
mu2e::GenParticlePrinter::Print(const mu2e::GenParticle& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;

  os << std::setiosflags(std::ios::fixed | std::ios::right);
  if(ind>=0) os << std::setw(4) << ind;

  os 
    << " " << std::setw(5) << obj.pdgId()
    << "  "
    << " " << std::setw(8) << std::setprecision(1) << obj.position().x()
    << " " << std::setw(8) << std::setprecision(1) << obj.position().y()
    << " " << std::setw(8) << std::setprecision(1) << obj.position().z()
    << "  "
    << " " << std::setw(8) << std::setprecision(1) << obj.momentum().x()
    << " " << std::setw(8) << std::setprecision(1) << obj.momentum().y()
    << " " << std::setw(8) << std::setprecision(1) << obj.momentum().z()
    << " " << std::setw(7) << std::setprecision(1) << obj.time()
    << " " << std::setw(7) << std::setprecision(1) << obj.properTime()
    << "   "
    << " " << std::setiosflags(std::ios::left) << obj.generatorId().name()
    << std::endl;

}

void 
mu2e::GenParticlePrinter::PrintHeader(const std::string& tag, std::ostream& os) {
  if(verbose()<1) return;
  os << "\nProductPrint " << tag << "\n";
}

void 
mu2e::GenParticlePrinter::PrintListHeader(std::ostream& os) {
  if(verbose()<1) return;
  os << "ind   pdgId            Position                     Momentum            time   ptime            name\n";

}

