
#include "Print/inc/BkgQualPrinter.hh"
#include "art/Framework/Principal/Provenance.h"
#include <string>
#include <iomanip>

void 
mu2e::BkgQualPrinter::Print(art::Event const& event,
				std::ostream& os) {
  if(verbose()<1) return;
  if(tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector< art::Handle<BkgQualCollection> > vah = event.getMany<BkgQualCollection>();
    for (auto const & ah : vah) Print(ah);
  } else {
    // print requested instances
    for(const auto& tag : tags() ) {
      auto ih = event.getValidHandle<BkgQualCollection>(tag);
      Print(ih);
    }
  }
}

void 
mu2e::BkgQualPrinter::Print(const art::Handle<BkgQualCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::BkgQualPrinter::Print(const art::ValidHandle<BkgQualCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::BkgQualPrinter::Print(const BkgQualCollection& coll, std::ostream& os) {
  if(verbose()<1) return;
  os << "BkgQualCollection has " << coll.size() << " clusters\n";
  if(!coll.empty()) {
    auto const& x = coll[0];
    os << "variables: " << std::endl;
    for(auto i: x.varNames()) {
      os << "     " << std::setw(3) << i.second 
	 << "   " << i.first << std::endl;
    }
  }
  if(verbose()==1) PrintListHeader();
  int i = 0;
  for(const auto& obj: coll) Print(obj, i++);
}

void 
mu2e::BkgQualPrinter::Print(const art::Ptr<BkgQual>& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;
  Print(*obj,ind);
}

void 
mu2e::BkgQualPrinter::Print(const mu2e::BkgQual& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;

  os << std::setiosflags(std::ios::fixed | std::ios::right);
  if(ind>=0) os << std::setw(4) << ind;

  if(verbose()==1) {
    os 
      << " " << std::setw(8) << std::setprecision(3) << obj.MVAOutput()
      << " " << std::setw(5) << (int)obj.status()
      << std::endl;
  }  

}

void 
mu2e::BkgQualPrinter::PrintHeader(const std::string& tag, std::ostream& os) {
  if(verbose()<1) return;
  os << "\nProductPrint " << tag << "\n";
}

void 
mu2e::BkgQualPrinter::PrintListHeader(std::ostream& os) {
  if(verbose()<1) return;
  os << "ind     mva   status \n";

}

