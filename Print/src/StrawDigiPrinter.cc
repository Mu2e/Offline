
#include "Print/inc/StrawDigiPrinter.hh"
#include "art/Framework/Principal/Provenance.h"
#include <string>
#include <iomanip>

void 
mu2e::StrawDigiPrinter::Print(art::Event const& event,
				std::ostream& os) {
  if(verbose()<1) return;
  if(tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector< art::Handle<StrawDigiCollection> > vah;
    event.getManyByType(vah);
    for (auto const & ah : vah) Print(ah);
  } else {
    // print requested instances
    for(const auto& tag : tags() ) {
      auto ih = event.getValidHandle<StrawDigiCollection>(tag);
      Print(ih);
    }
  }
}

void 
mu2e::StrawDigiPrinter::Print(const art::Handle<StrawDigiCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::StrawDigiPrinter::Print(const art::ValidHandle<StrawDigiCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::StrawDigiPrinter::Print(const StrawDigiCollection& coll, std::ostream& os) {
  if(verbose()<1) return;
  os << "StrawDigiCollection has " << coll.size() << " hits\n";
  if(verbose()==1) PrintListHeader();
  int i = 0;
  for(const auto& obj: coll) Print(obj, i++);
}

void 
mu2e::StrawDigiPrinter::Print(const art::Ptr<StrawDigi>& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;
  Print(*obj,ind);
}

void 
mu2e::StrawDigiPrinter::Print(const mu2e::StrawDigi& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;

  os << std::setiosflags(std::ios::fixed | std::ios::right);
  if(ind>=0) os << std::setw(4) << ind;

  os 
    << " " << std::setw(5) << obj.strawId().asUint16()
    << " " 
    << " " << std::setw(6) << obj.TDC()[0]
    << " " << std::setw(6) << obj.TDC()[1]
    << " " << std::setw(6) << obj.adcWaveform().size()
    << " " ;
  for(auto& i: obj.adcWaveform()) {
    os << " " << std::setw(4) << i;
  }
  os << std::endl;

}

void 
mu2e::StrawDigiPrinter::PrintHeader(const std::string& tag, std::ostream& os) {
  if(verbose()<1) return;
  os << "\nProductPrint " << tag << "\n";
}

void 
mu2e::StrawDigiPrinter::PrintListHeader(std::ostream& os) {
  if(verbose()<1) return;
  os << " ind StrwInd TDC0   TDC1    NADC  ADC\n";

}

