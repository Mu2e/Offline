
#include "Print/inc/CrvDigiMCPrinter.hh"
#include "art/Framework/Principal/Provenance.h"
#include <string>
#include <iomanip>

void 
mu2e::CrvDigiMCPrinter::Print(art::Event const& event,
				std::ostream& os) {
  if(verbose()<1) return;
  if(tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector< art::Handle<CrvDigiMCCollection> > vah = event.getMany<CrvDigiMCCollection>();
    for (auto const & ah : vah) Print(ah);
  } else {
    // print requested instances
    for(const auto& tag : tags() ) {
      auto ih = event.getValidHandle<CrvDigiMCCollection>(tag);
      Print(ih);
    }
  }
}

void 
mu2e::CrvDigiMCPrinter::Print(const art::Handle<CrvDigiMCCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::CrvDigiMCPrinter::Print(const art::ValidHandle<CrvDigiMCCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::CrvDigiMCPrinter::Print(const CrvDigiMCCollection& coll, std::ostream& os) {
  if(verbose()<1) return;
  os << "CrvDigiMCCollection has " << coll.size() << " hits\n";
  if(verbose()==1) PrintListHeader();
  int i = 0;
  for(const auto& obj: coll) Print(obj, i++);
}

void 
mu2e::CrvDigiMCPrinter::Print(const art::Ptr<CrvDigiMC>& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;
  Print(*obj,ind);
}

void 
mu2e::CrvDigiMCPrinter::Print(const mu2e::CrvDigiMC& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;

  os << std::setiosflags(std::ios::fixed | std::ios::right);
  if(ind>=0) os << std::setw(4) << ind;

  int simp = -1;
  if(obj.GetSimParticle()) {
    simp = int(obj.GetSimParticle().key());
  }

  os 
    << " " << std::setw(5) << obj.GetScintillatorBarIndex()
    << " " << std::setw(5) << obj.GetSiPMNumber()
    << " " << std::setw(5) << simp 
    << " " << std::setw(5) << obj.GetCrvSteps().size()
    << "   " << std::setw(8) << std::setprecision(3) << obj.GetStartTime();
  os << std::endl;

}

void 
mu2e::CrvDigiMCPrinter::PrintHeader(const std::string& tag, std::ostream& os) {
  if(verbose()<1) return;
  os << "\nProductPrint " << tag << "\n";
}

void 
mu2e::CrvDigiMCPrinter::PrintListHeader(std::ostream& os) {
  if(verbose()<1) return;
  os << "ind    Bar  SiPM   SimP   NStep   time\n";


}

