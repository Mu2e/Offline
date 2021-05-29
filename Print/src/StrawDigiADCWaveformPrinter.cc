
#include "Print/inc/StrawDigiADCWaveformPrinter.hh"
#include "art/Framework/Principal/Provenance.h"
#include <string>
#include <iomanip>

void 
mu2e::StrawDigiADCWaveformPrinter::Print(art::Event const& event,
				std::ostream& os) {
  if(verbose()<1) return;
  if(tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector< art::Handle<StrawDigiADCWaveformCollection> > vah = event.getMany<StrawDigiADCWaveformCollection>();
    for (auto const & ah : vah) Print(ah);
  } else {
    // print requested instances
    for(const auto& tag : tags() ) {
      auto ih = event.getValidHandle<StrawDigiADCWaveformCollection>(tag);
      Print(ih);
    }
  }
}

void 
mu2e::StrawDigiADCWaveformPrinter::Print(const art::Handle<StrawDigiADCWaveformCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::StrawDigiADCWaveformPrinter::Print(const art::ValidHandle<StrawDigiADCWaveformCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::StrawDigiADCWaveformPrinter::Print(const StrawDigiADCWaveformCollection& coll, std::ostream& os) {
  if(verbose()<1) return;
  os << "StrawDigiADCWaveformCollection has " << coll.size() << " hits\n";
  if(verbose()==1) PrintListHeader();
  int i = 0;
  for(const auto& obj: coll) Print(obj, i++);
}

void 
mu2e::StrawDigiADCWaveformPrinter::Print(const art::Ptr<StrawDigiADCWaveform>& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;
  Print(*obj,ind);
}

void 
mu2e::StrawDigiADCWaveformPrinter::Print(const mu2e::StrawDigiADCWaveform& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;

  os << std::setiosflags(std::ios::fixed | std::ios::right);
  if(ind>=0) os << std::setw(4) << ind;

  for (auto &i : obj.samples()){
    os << " " << std::setw(6) << i;
  }
  os << std::endl;

}

void 
mu2e::StrawDigiADCWaveformPrinter::PrintHeader(const std::string& tag, std::ostream& os) {
  if(verbose()<1) return;
  os << "\nProductPrint " << tag << "\n";
}

void 
mu2e::StrawDigiADCWaveformPrinter::PrintListHeader(std::ostream& os) {
  if(verbose()<1) return;
  os << " ind Waveform\n";

}

