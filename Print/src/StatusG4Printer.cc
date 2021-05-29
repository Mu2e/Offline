
#include "Print/inc/StatusG4Printer.hh"
#include "art/Framework/Principal/Provenance.h"
#include <string>
#include <iomanip>

void 
mu2e::StatusG4Printer::Print(art::Event const& event,
				std::ostream& os) {
  if(verbose()<1) return;
  if(tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector< art::Handle<StatusG4> > vah = event.getMany<StatusG4>();
    for (auto const & ah : vah) Print(ah);
  } else {
    // print requested instances
    for(const auto& tag : tags() ) {
      auto ih = event.getValidHandle<StatusG4>(tag);
      Print(ih);
    }
  }
}

void 
mu2e::StatusG4Printer::Print(const art::Handle<StatusG4>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::StatusG4Printer::Print(const art::ValidHandle<StatusG4>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}


void 
mu2e::StatusG4Printer::Print(const mu2e::StatusG4& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;

  os << std::setiosflags(std::ios::fixed | std::ios::right);

  os 
    << "  status: " << std::setw(3) << obj.status()
    << " nG4Tracks: " << std::setw(7) << obj.nG4Tracks()
    << " overflowSimP: " << std::setw(2) << obj.overflowSimParticles()
    << " nStepLimit: " << std::setw(2) << obj.nKilledStepLimit()
    << " CPU: " << std::setw(8) << std::setprecision(3) << obj.cpuTime()
    << " Wall: " << std::setw(8) << std::setprecision(3) <<  obj.realTime();
  os << std::endl;

}

void 
mu2e::StatusG4Printer::PrintHeader(const std::string& tag, std::ostream& os) {
  if(verbose()<1) return;
  os << "\nProductPrint " << tag << "\n";
}


