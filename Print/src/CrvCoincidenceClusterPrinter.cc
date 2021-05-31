
#include "Print/inc/CrvCoincidenceClusterPrinter.hh"
#include "art/Framework/Principal/Provenance.h"
#include <string>
#include <iomanip>

void 
mu2e::CrvCoincidenceClusterPrinter::Print(art::Event const& event,
				std::ostream& os) {
  if(verbose()<1) return;
  if(tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector< art::Handle<CrvCoincidenceClusterCollection> > vah = event.getMany<CrvCoincidenceClusterCollection>();
    for (auto const & ah : vah) Print(ah);
  } else {
    // print requested instances
    for(const auto& tag : tags() ) {
      auto ih = event.getValidHandle<CrvCoincidenceClusterCollection>(tag);
      Print(ih);
    }
  }
}

void 
mu2e::CrvCoincidenceClusterPrinter::Print(const art::Handle<CrvCoincidenceClusterCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::CrvCoincidenceClusterPrinter::Print(const art::ValidHandle<CrvCoincidenceClusterCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::CrvCoincidenceClusterPrinter::Print(const CrvCoincidenceClusterCollection& coll, std::ostream& os) {
  if(verbose()<1) return;
  os << "CrvCoincidenceClusterCollection has " << coll.size() << " clusters\n";
  if(verbose()==1) PrintListHeader();
  int i = 0;
  for(const auto& obj: coll) Print(obj, i++);
}

void 
mu2e::CrvCoincidenceClusterPrinter::Print(const art::Ptr<CrvCoincidenceCluster>& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;
  Print(*obj,ind);
}

void 
mu2e::CrvCoincidenceClusterPrinter::Print(const mu2e::CrvCoincidenceCluster& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;

  os << std::setiosflags(std::ios::fixed | std::ios::right);
  if(ind>=0) os << std::setw(4) << ind;

  os 
    << " " << std::setw(5) << obj.GetCrvSectorType()
    << " " << std::setw(8) << std::setprecision(1) << obj.GetPEs()
    << " " << std::setw(8) << std::setprecision(1) << obj.GetStartTime()
    << " " << std::setw(8) << std::setprecision(1) << obj.GetEndTime()
    << "   " 
    << " " << std::setw(8) << std::setprecision(1) << obj.GetAvgCounterPos().x()
    << " " << std::setw(8) << std::setprecision(1) << obj.GetAvgCounterPos().y()
    << " " << std::setw(8) << std::setprecision(1) << obj.GetAvgCounterPos().z()
    << " ";
  os << std::endl;
  if(verbose()>1) {
    os << "    RecoPulses:";
    for(auto pp : obj.GetCrvRecoPulses()) {
      os << " " << pp.key();
    }
    os << std::endl;
  }

}

void 
mu2e::CrvCoincidenceClusterPrinter::PrintHeader(const std::string& tag, std::ostream& os) {
  if(verbose()<1) return;
  os << "\nProductPrint " << tag << "\n";
}

void 
mu2e::CrvCoincidenceClusterPrinter::PrintListHeader(std::ostream& os) {
  if(verbose()<1) return;
  os << "ind   SecType   PEs   t_start   t_end          avg position\n";

}

