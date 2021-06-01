
#include "Print/inc/MCTrajectoryPrinter.hh"
#include "art/Framework/Principal/Provenance.h"
#include <string>
#include <iomanip>

void 
mu2e::MCTrajectoryPrinter::Print(art::Event const& event,
				std::ostream& os) {
  if(verbose()<1) return;
  if(tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector< art::Handle<MCTrajectoryCollection> > vah = event.getMany<MCTrajectoryCollection>();
    for (auto const & ah : vah) Print(ah);
  } else {
    // print requested instances
    for(const auto& tag : tags() ) {
      auto ih = event.getValidHandle<MCTrajectoryCollection>(tag);
      Print(ih);
    }
  }
}

void 
mu2e::MCTrajectoryPrinter::Print(const art::Handle<MCTrajectoryCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::MCTrajectoryPrinter::Print(const art::ValidHandle<MCTrajectoryCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::MCTrajectoryPrinter::Print(const MCTrajectoryCollection& coll, std::ostream& os) {
  if(verbose()<1) return;
  os << "MCTrajectoryCollection has " << coll.size() << " trajectories\n";
  if(verbose()==1) PrintListHeader();
  int i = 0;
  for(const auto& obj: coll) Print(obj.second, i++);
}

void 
mu2e::MCTrajectoryPrinter::Print(const art::Ptr<MCTrajectory>& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;
  Print(*obj,ind);
}

void 
mu2e::MCTrajectoryPrinter::Print(const mu2e::MCTrajectory& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;

  int pkey = obj.simid();
  const auto& points = obj.points();

  os << std::setiosflags(std::ios::fixed | std::ios::right);
  if(ind>=0) os << std::setw(4) << ind;

  if(verbose()==1) {


    const auto& v0 = points.front();
    const auto& v1 = points.back();

    os 
      << " " << std::setw(8) << pkey
      << " " << std::setw(4) << points.size()
      << "  "
      << " " << std::setw(8) << std::setprecision(1) << v0.x()
      << " " << std::setw(8) << std::setprecision(1) << v0.y()
      << " " << std::setw(8) << std::setprecision(1) << v0.z()
      << " " << std::setw(9) << std::setprecision(1) << v0.t()
      << " " << std::setw(9) << std::setprecision(1) << v0.kineticEnergy()
      << "  "
      << " " << std::setw(8) << std::setprecision(1) << v1.x()
      << " " << std::setw(8) << std::setprecision(1) << v1.y()
      << " " << std::setw(8) << std::setprecision(1) << v1.z()
      << " " << std::setw(9) << std::setprecision(1) << v1.t()
      << " " << std::setw(9) << std::setprecision(1) << v1.kineticEnergy()
      << std::endl;

  } else {

    os 
      << " parentKey: " << std::setw(8) << pkey
      << "  npoint: " << std::setw(4) << points.size() << "\n";
    int i=0;
    for(auto const& pp : points ) {
      os
	<< "  " << std::setw(4) << i++ 
	<< " " << std::setw(8) << std::setprecision(1) << pp.x()
	<< " " << std::setw(8) << std::setprecision(1) << pp.y()
	<< " " << std::setw(8) << std::setprecision(1) << pp.z()
	<< " " << std::setw(9) << std::setprecision(1) << pp.t()
	<< " " << std::setw(9) << std::setprecision(1) << pp.kineticEnergy()
        << std::endl;
    }

  } // end if verbose

}

void 
mu2e::MCTrajectoryPrinter::PrintHeader(const std::string& tag, std::ostream& os) {
  if(verbose()<1) return;
  os << "\nProductPrint " << tag << "\n";
}

void 
mu2e::MCTrajectoryPrinter::PrintListHeader(std::ostream& os) {
  if(verbose()<1) return;
  os << "ind   parent  npoint       first Point            firstT     firstEk            last Point           lastT    lastEk\n";
}

