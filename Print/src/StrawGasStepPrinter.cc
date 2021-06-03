
#include "Print/inc/StrawGasStepPrinter.hh"
#include "art/Framework/Principal/Provenance.h"
#include <string>
#include <iomanip>

void 
mu2e::StrawGasStepPrinter::Print(art::Event const& event,
				std::ostream& os) {
  if(verbose()<1) return;
  if(tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector< art::Handle<StrawGasStepCollection> > vah = event.getMany<StrawGasStepCollection>();
    for (auto const & ah : vah) Print(ah);
  } else {
    // print requested instances
    for(const auto& tag : tags() ) {
      auto ih = event.getValidHandle<StrawGasStepCollection>(tag);
      Print(ih);
    }
  }
}

void 
mu2e::StrawGasStepPrinter::Print(const art::Handle<StrawGasStepCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::StrawGasStepPrinter::Print(const art::ValidHandle<StrawGasStepCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::StrawGasStepPrinter::Print(const StrawGasStepCollection& coll, std::ostream& os) {
  if(verbose()<1) return;
  os << "StrawGasStepCollection has " << coll.size() << " steps\n";
  if(verbose()==1) PrintListHeader();
  int i = 0;
  for(const auto& obj: coll) Print(obj, i++);
}

void 
mu2e::StrawGasStepPrinter::Print(const art::Ptr<StrawGasStep>& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;
  Print(*obj,ind);
}

void 
mu2e::StrawGasStepPrinter::Print(const mu2e::StrawGasStep& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;

  os << std::setiosflags(std::ios::fixed | std::ios::right);
  if(ind>=0) os << std::setw(4) << ind;

  long unsigned int pkey=0;
  auto const& pptr = obj.simParticle();
  if(pptr) pkey = pptr->id().asUint();

  os 
    << " " << std::setw(5) << pkey
    << " " << std::setw(8) << obj.strawId().asUint16()
    << " " << std::setw(10) << std::setprecision(6) << obj.ionizingEdep()
    << " " << std::setw(8) << std::setprecision(3) << obj.stepLength()
    << " " << std::setw(8) << std::setprecision(1) << obj.time()
    << std::endl;

}

void 
mu2e::StrawGasStepPrinter::PrintHeader(const std::string& tag, std::ostream& os) {
  if(verbose()<1) return;
  os << "\nProductPrint " << tag << "\n";
}

void 
mu2e::StrawGasStepPrinter::PrintListHeader(std::ostream& os) {
  if(verbose()<1) return;
  os << " ind SimPart StrwInd    eDep     length    time\n";

}

