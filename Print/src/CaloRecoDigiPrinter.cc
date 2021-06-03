
#include "Print/inc/CaloRecoDigiPrinter.hh"
#include "art/Framework/Principal/Provenance.h"
#include <string>
#include <iomanip>

void 
mu2e::CaloRecoDigiPrinter::Print(art::Event const& event,
				std::ostream& os) {
  if(verbose()<1) return;
  if(tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector< art::Handle<CaloRecoDigiCollection> > vah = event.getMany<CaloRecoDigiCollection>();
    for (auto const & ah : vah) Print(ah);
  } else {
    // print requested instances
    for(const auto& tag : tags() ) {
      auto ih = event.getValidHandle<CaloRecoDigiCollection>(tag);
      Print(ih);
    }
  }
}

void 
mu2e::CaloRecoDigiPrinter::Print(const art::Handle<CaloRecoDigiCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::CaloRecoDigiPrinter::Print(const art::ValidHandle<CaloRecoDigiCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::CaloRecoDigiPrinter::Print(const CaloRecoDigiCollection& coll, std::ostream& os) {
  if(verbose()<1) return;
  os << "CaloRecoDigiCollection has " << coll.size() << " hits\n";
  if(verbose()==1) PrintListHeader();
  int i = 0;
  for(const auto& obj: coll) Print(obj, i++);
}

void 
mu2e::CaloRecoDigiPrinter::Print(const art::Ptr<CaloRecoDigi>& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;
  Print(*obj,ind);
}

void 
mu2e::CaloRecoDigiPrinter::Print(const mu2e::CaloRecoDigi& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;

  os << std::setiosflags(std::ios::fixed | std::ios::right);
  if(ind>=0) os << std::setw(4) << ind;

  os 
    << " " << std::setw(5) << obj.SiPMID()
    << " " 
    << " " << std::setw(8) << std::setprecision(1) << obj.energyDep()
    << " " << std::setw(8) << std::setprecision(1) << obj.energyDepErr()
    << " " << std::setw(8) << std::setprecision(2) << obj.time()
    << " " << std::setw(8) << std::setprecision(2) << obj.timeErr()
    << " " << std::setw(8) << std::setprecision(1) << obj.chi2()
    << " " << std::setw(5) << obj.ndf()
    << " " << std::setw(3) << obj.pileUp()
    << std::endl;
 
}

void 
mu2e::CaloRecoDigiPrinter::PrintHeader(const std::string& tag, std::ostream& os) {
  if(verbose()<1) return;
  os << "\nProductPrint " << tag << "\n";
}

void 
mu2e::CaloRecoDigiPrinter::PrintListHeader(std::ostream& os) {
  if(verbose()<1) return;
  os << "ind   SiPMID     energy   e_err    time     t_err     chi2   ndf  pileup\n";

}


