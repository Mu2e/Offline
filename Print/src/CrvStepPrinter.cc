
#include "Print/inc/CrvStepPrinter.hh"
#include "art/Framework/Principal/Provenance.h"
#include <string>
#include <iomanip>

void 
mu2e::CrvStepPrinter::Print(art::Event const& event,
				std::ostream& os) {
  if(verbose()<1) return;
  if(tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector< art::Handle<CrvStepCollection> > vah = event.getMany<CrvStepCollection>();
    for (auto const & ah : vah) Print(ah);
  } else {
    // print requested instances
    for(const auto& tag : tags() ) {
      auto ih = event.getValidHandle<CrvStepCollection>(tag);
      Print(ih);
    }
  }
}

void 
mu2e::CrvStepPrinter::Print(const art::Handle<CrvStepCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::CrvStepPrinter::Print(const art::ValidHandle<CrvStepCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::CrvStepPrinter::Print(const CrvStepCollection& coll, std::ostream& os) {
  if(verbose()<1) return;
  os << "CrvStepCollection has " << coll.size() << " steps\n";
  if(verbose()==1) PrintListHeader();
  int i = 0;
  for(const auto& obj: coll) Print(obj, i++);
}

void 
mu2e::CrvStepPrinter::Print(const art::Ptr<CrvStep>& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;
  Print(*obj,ind);
}

void 
mu2e::CrvStepPrinter::Print(const mu2e::CrvStep& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;

  if( obj.visibleEDep() < _eCut ) return;

  os << std::setiosflags(std::ios::fixed | std::ios::right);
  if(ind>=0) os << std::setw(4) << ind;

  auto const& ptr = obj.simParticle();
  int pkey = -1;
  if(ptr) pkey = ptr.key();

  auto const & pos = obj.startPosition();
  if(verbose()==1) {
    os 
      << " " << std::setw(5) << obj.barIndex()
      << " " << std::setw(9) << pkey
      << " " << std::setw(6) << std::setprecision(1) << obj.visibleEDep()
      << " " << std::setw(8) << std::setprecision(1) << obj.startTime()
      << " " << std::setw(9) << std::setprecision(2) << pos.x()
      << " " << std::setw(9) << std::setprecision(2) << pos.y()
      << " " << std::setw(9) << std::setprecision(2) << pos.z()
      << " " << std::setw(8) << std::setprecision(1) << obj.startMomentum().mag() 
      << "  " 
      << std::endl;
  }
}

void 
mu2e::CrvStepPrinter::PrintHeader(const std::string& tag, std::ostream& os) {
  if(verbose()<1) return;
  os << "\nProductPrint " << tag << "\n";
}

void 
mu2e::CrvStepPrinter::PrintListHeader(std::ostream& os) {
  if(verbose()<1) return;
  os << "ind    bar    simP     EDep   startT           startPos              startP\n";

}

