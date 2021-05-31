
#include "Print/inc/CaloDigiPrinter.hh"
#include "art/Framework/Principal/Provenance.h"
#include <string>
#include <iomanip>

void 
mu2e::CaloDigiPrinter::Print(art::Event const& event,
				std::ostream& os) {
  if(verbose()<1) return;
  if(tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector< art::Handle<CaloDigiCollection> > vah = event.getMany<CaloDigiCollection>();
    for (auto const & ah : vah) Print(ah);
  } else {
    // print requested instances
    for(const auto& tag : tags() ) {
      auto ih = event.getValidHandle<CaloDigiCollection>(tag);
      Print(ih);
    }
  }
}

void 
mu2e::CaloDigiPrinter::Print(const art::Handle<CaloDigiCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::CaloDigiPrinter::Print(const art::ValidHandle<CaloDigiCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::CaloDigiPrinter::Print(const CaloDigiCollection& coll, std::ostream& os) {
  if(verbose()<1) return;
  os << "CaloDigiCollection has " << coll.size() << " hits\n";
  if(verbose()==1) PrintListHeader();
  int i = 0;
  for(const auto& obj: coll) Print(obj, i++);
}

void 
mu2e::CaloDigiPrinter::Print(const art::Ptr<CaloDigi>& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;
  Print(*obj,ind);
}

void 
mu2e::CaloDigiPrinter::Print(const mu2e::CaloDigi& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;

  os << std::setiosflags(std::ios::fixed | std::ios::right);
  if(ind>=0) os << std::setw(4) << ind;

  os 
    << " " << std::setw(5) << obj.SiPMID()
    << " " << std::setw(5) << obj.t0()
    << " ";
  for(auto i: obj.waveform()) {
    os << " " << std::setw(5) << i;
  }
  os << std::endl;

}

void 
mu2e::CaloDigiPrinter::PrintHeader(const std::string& tag, std::ostream& os) {
  if(verbose()<1) return;
  os << "\nProductPrint " << tag << "\n";
}

void 
mu2e::CaloDigiPrinter::PrintListHeader(std::ostream& os) {
  if(verbose()<1) return;
  os << "ind   roId   time     waveform\n";

}

