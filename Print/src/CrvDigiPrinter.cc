
#include "Print/inc/CrvDigiPrinter.hh"
#include "art/Framework/Principal/Provenance.h"
#include <string>
#include <iomanip>

void 
mu2e::CrvDigiPrinter::Print(art::Event const& event,
				std::ostream& os) {
  if(verbose()<1) return;
  if(tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector< art::Handle<CrvDigiCollection> > vah;
    event.getManyByType(vah);
    for (auto const & ah : vah) Print(ah);
  } else {
    // print requested instances
    for(const auto& tag : tags() ) {
      auto ih = event.getValidHandle<CrvDigiCollection>(tag);
      Print(ih);
    }
  }
}

void 
mu2e::CrvDigiPrinter::Print(const art::Handle<CrvDigiCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::CrvDigiPrinter::Print(const art::ValidHandle<CrvDigiCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::CrvDigiPrinter::Print(const CrvDigiCollection& coll, std::ostream& os) {
  if(verbose()<1) return;
  os << "CrvDigiCollection has " << coll.size() << " hits\n";
  if(verbose()==1) PrintListHeader();
  int i = 0;
  for(const auto& obj: coll) Print(obj, i++);
}

void 
mu2e::CrvDigiPrinter::Print(const art::Ptr<CrvDigi>& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;
  Print(*obj,ind);
}

void 
mu2e::CrvDigiPrinter::Print(const mu2e::CrvDigi& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;

  os << std::setiosflags(std::ios::fixed | std::ios::right);
  if(ind>=0) os << std::setw(4) << ind;

  os 
    << " " << std::setw(5) << obj.GetScintillatorBarIndex()
    << " " << std::setw(5) << obj.GetSiPMNumber()
    << " " << std::setw(5) << obj.GetStartTDC()
    << " ";
  for(auto i: obj.GetADCs()) {
    os << " " << std::setw(5) << i;
  }
  os << std::endl;

}

void 
mu2e::CrvDigiPrinter::PrintHeader(const std::string& tag, std::ostream& os) {
  if(verbose()<1) return;
  os << "\nProductPrint " << tag << "\n";
}

void 
mu2e::CrvDigiPrinter::PrintListHeader(std::ostream& os) {
  if(verbose()<1) return;
  os << "ind    Bar  SiPM   TDC     ADC waveform\n";

}


