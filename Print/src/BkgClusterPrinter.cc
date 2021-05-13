
#include "Print/inc/BkgClusterPrinter.hh"
#include "art/Framework/Principal/Provenance.h"
#include <string>
#include <iomanip>

void 
mu2e::BkgClusterPrinter::Print(art::Event const& event,
				std::ostream& os) {
  if(verbose()<1) return;
  if(tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector< art::Handle<BkgClusterCollection> > vah = event.getMany<BkgClusterCollection>();
    for (auto const & ah : vah) Print(ah);
  } else {
    // print requested instances
    for(const auto& tag : tags() ) {
      auto ih = event.getValidHandle<BkgClusterCollection>(tag);
      Print(ih);
    }
  }
}

void 
mu2e::BkgClusterPrinter::Print(const art::Handle<BkgClusterCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::BkgClusterPrinter::Print(const art::ValidHandle<BkgClusterCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::BkgClusterPrinter::Print(const BkgClusterCollection& coll, std::ostream& os) {
  if(verbose()<1) return;
  os << "BkgClusterCollection has " << coll.size() << " clusters\n";
  if(verbose()==1) PrintListHeader();
  int i = 0;
  for(const auto& obj: coll) Print(obj, i++);
}

void 
mu2e::BkgClusterPrinter::Print(const art::Ptr<BkgCluster>& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;
  Print(*obj,ind);
}

void 
mu2e::BkgClusterPrinter::Print(const mu2e::BkgCluster& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;

  os << std::setiosflags(std::ios::fixed | std::ios::right);
  if(ind>=0) os << std::setw(4) << ind;

  if(verbose()==1) {
    os 
      << " " << std::setw(8) << std::setprecision(3) << obj.pos().x()
      << " " << std::setw(8) << std::setprecision(3) << obj.pos().y()
      << " " << std::setw(8) << std::setprecision(3) << obj.pos().z()
      << "   " << std::setw(7) << std::setprecision(1) << obj.time()
      << " " << std::setw(5) << obj.hits().size()
      << std::endl;
  } else if(verbose()==2) {
    os 
      << "   pos: " << std::setw(9) << std::setprecision(3) << obj.pos().x()
      << " " << std::setw(9) << std::setprecision(3) << obj.pos().y()
      << " " << std::setw(8) << std::setprecision(3) << obj.pos().z()
      << "   time: " << std::setw(7) << std::setprecision(1) << obj.time()

      << "  hits: " << std::setw(4) << obj.hits().size()
      << "\n";

    os 
      << "   flag: " ;
    for(auto sn: obj.flag().bitNames()) { 
      if(obj.flag().hasAnyProperty(BkgClusterFlag(sn.first))) 
	os << " " << sn.first;
    }
    os 
      << "\n";

    os 
      << "   hits: " ;
    for(auto bch : obj.hits()) {
      os << " " << std::setw(5) << bch;
    }
    os << std::endl;
  }
}

void 
mu2e::BkgClusterPrinter::PrintHeader(const std::string& tag, std::ostream& os) {
  if(verbose()<1) return;
  os << "\nProductPrint " << tag << "\n";
}

void 
mu2e::BkgClusterPrinter::PrintListHeader(std::ostream& os) {
  if(verbose()<1) return;
  os << "ind           position             time    nhits\n";

}

