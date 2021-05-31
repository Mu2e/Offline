
#include "Print/inc/CrvRecoPulsePrinter.hh"
#include "art/Framework/Principal/Provenance.h"
#include <string>
#include <iomanip>

void 
mu2e::CrvRecoPulsePrinter::Print(art::Event const& event,
				std::ostream& os) {
  if(verbose()<1) return;
  if(tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector< art::Handle<CrvRecoPulseCollection> > vah = event.getMany<CrvRecoPulseCollection>();
    for (auto const & ah : vah) Print(ah);
  } else {
    // print requested instances
    for(const auto& tag : tags() ) {
      auto ih = event.getValidHandle<CrvRecoPulseCollection>(tag);
      Print(ih);
    }
  }
}

void 
mu2e::CrvRecoPulsePrinter::Print(const art::Handle<CrvRecoPulseCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::CrvRecoPulsePrinter::Print(const art::ValidHandle<CrvRecoPulseCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::CrvRecoPulsePrinter::Print(const CrvRecoPulseCollection& coll, std::ostream& os) {
  if(verbose()<1) return;
  os << "CrvRecoPulseCollection has " << coll.size() << " hits\n";
  if(verbose()==1) PrintListHeader();
  int i = 0;
  for(const auto& obj: coll) Print(obj, i++);
}

void 
mu2e::CrvRecoPulsePrinter::Print(const art::Ptr<CrvRecoPulse>& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;
  Print(*obj,ind);
}

void 
mu2e::CrvRecoPulsePrinter::Print(const mu2e::CrvRecoPulse& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;

  os << std::setiosflags(std::ios::fixed | std::ios::right);
  if(ind>=0) os << std::setw(4) << ind;

  os 
    << " " << std::setw(5) << obj.GetScintillatorBarIndex()
    << " " << std::setw(5) << obj.GetSiPMNumber()
    << " " << std::setw(5) << obj.GetPEs()
    << " " << std::setw(5) << obj.GetPEsPulseHeight()
    << " " << std::setw(10) << std::setprecision(1) << obj.GetPulseTime()
    << " " << std::setw(10) << std::setprecision(3) << obj.GetPulseFitChi2()
    << " " << std::setw(9) << std::setprecision(1) << obj.GetLEtime();
  os << std::endl;

}

void 
mu2e::CrvRecoPulsePrinter::PrintHeader(const std::string& tag, std::ostream& os) {
  if(verbose()<1) return;
  os << "\nProductPrint " << tag << "\n";
}

void 
mu2e::CrvRecoPulsePrinter::PrintListHeader(std::ostream& os) {
  if(verbose()<1) return;
  os << "ind    Bar  SiPM   PEs PEheight    time     chi2     LEtime\n";

}

