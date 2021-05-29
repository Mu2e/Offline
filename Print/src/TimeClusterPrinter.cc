
#include "Print/inc/TimeClusterPrinter.hh"
#include "art/Framework/Principal/Provenance.h"
#include <string>
#include <iomanip>

void 
mu2e::TimeClusterPrinter::Print(art::Event const& event,
				std::ostream& os) {
  if(verbose()<1) return;
  if(tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector< art::Handle<TimeClusterCollection> > vah = event.getMany<TimeClusterCollection>();
    for (auto const & ah : vah) Print(ah);
  } else {
    // print requested instances
    for(const auto& tag : tags() ) {
      auto ih = event.getValidHandle<TimeClusterCollection>(tag);
      Print(ih);
    }
  }
}

void 
mu2e::TimeClusterPrinter::Print(const art::Handle<TimeClusterCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::TimeClusterPrinter::Print(const art::ValidHandle<TimeClusterCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::TimeClusterPrinter::Print(const TimeClusterCollection& coll, std::ostream& os) {
  if(verbose()<1) return;
  os << "TimeClusterCollection has " << coll.size() << " hits\n";
  if(verbose()==1) PrintListHeader();
  int i = 0;
  for(const auto& obj: coll) Print(obj, i++);
}

void 
mu2e::TimeClusterPrinter::Print(const art::Ptr<TimeCluster>& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;
  Print(*obj,ind);
}

void 
mu2e::TimeClusterPrinter::Print(const mu2e::TimeCluster& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;

  os << std::setiosflags(std::ios::fixed | std::ios::right);
  if(ind>=0) os << std::setw(4) << ind;

  os 
    << " " << std::setw(5) << obj.nhits()
    << " " 
    << " " << std::setw(8) << std::setprecision(3) << obj.position().x()
    << " " << std::setw(8) << std::setprecision(3) << obj.position().y()
    << " " << std::setw(9) << std::setprecision(3) << obj.position().z()
    << " " << std::setw(7) << std::setprecision(1) << obj.t0().t0()
    << std::endl;

}

void 
mu2e::TimeClusterPrinter::PrintHeader(const std::string& tag, std::ostream& os) {
  if(verbose()<1) return;
  os << "\nProductPrint " << tag << "\n";
}

void 
mu2e::TimeClusterPrinter::PrintListHeader(std::ostream& os) {
  if(verbose()<1) return;
  os << " ind   nhits      x       y       z       time\n";

}

