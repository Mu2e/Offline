
#include "Print/inc/CaloShowerStepPrinter.hh"
#include "art/Framework/Principal/Provenance.h"
#include <string>
#include <iomanip>

void 
mu2e::CaloShowerStepPrinter::Print(art::Event const& event,
				std::ostream& os) {
  if(verbose()<1) return;
  if(tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector< art::Handle<CaloShowerStepCollection> > vah = event.getMany<CaloShowerStepCollection>();
    for (auto const & ah : vah) Print(ah);
  } else {
    // print requested instances
    for(const auto& tag : tags() ) {
      auto ih = event.getValidHandle<CaloShowerStepCollection>(tag);
      Print(ih);
    }
  }
}

void 
mu2e::CaloShowerStepPrinter::Print(const art::Handle<CaloShowerStepCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::CaloShowerStepPrinter::Print(const art::ValidHandle<CaloShowerStepCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::CaloShowerStepPrinter::Print(const CaloShowerStepCollection& coll, std::ostream& os) {
  if(verbose()<1) return;
  os << "CaloShowerStepCollection has " << coll.size() << " steps\n";
  if(verbose()==1) PrintListHeader();
  int i = 0;
  for(const auto& obj: coll) Print(obj, i++);
}

void 
mu2e::CaloShowerStepPrinter::Print(const art::Ptr<CaloShowerStep>& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;
  Print(*obj,ind);
}

void 
mu2e::CaloShowerStepPrinter::Print(const mu2e::CaloShowerStep& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;

  if( obj.energyDepG4() < _eCut ) return;

  os << std::setiosflags(std::ios::fixed | std::ios::right);
  if(ind>=0) os << std::setw(4) << ind;

  const CLHEP::Hep3Vector& pos = obj.position();
  if(verbose()==1) {
    os 
      << " " << std::setw(5) << obj.volumeG4ID()
      << " " << std::setw(12) << obj.simParticle().key()
      << " " << std::setw(8) << std::setprecision(2) << obj.energyDepG4()
      << " " << std::setw(8) << std::setprecision(1) << obj.time() 
      << "  " 
      << " " << std::setw(8) << std::setprecision(2) << pos.x()
      << " " << std::setw(8) << std::setprecision(2) << pos.y()
      << " " << std::setw(8) << std::setprecision(2) << pos.z()
      << std::endl;
  }
}

void 
mu2e::CaloShowerStepPrinter::PrintHeader(const std::string& tag, std::ostream& os) {
  if(verbose()<1) return;
  os << "\nProductPrint " << tag << "\n";
}

void 
mu2e::CaloShowerStepPrinter::PrintListHeader(std::ostream& os) {
  if(verbose()<1) return;
  os << "ind    vol        SimPart   energy    time            position\n";

}

