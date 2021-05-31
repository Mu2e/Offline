
#include "Print/inc/TrackClusterMatchPrinter.hh"
#include "art/Framework/Principal/Provenance.h"
#include <string>
#include <iomanip>

void 
mu2e::TrackClusterMatchPrinter::Print(art::Event const& event,
				std::ostream& os) {
  if(verbose()<1) return;
  if(tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector< art::Handle<TrackClusterMatchCollection> > vah = event.getMany<TrackClusterMatchCollection>();
    for (auto const & ah : vah) Print(ah);
  } else {
    // print requested instances
    for(const auto& tag : tags() ) {
      auto ih = event.getValidHandle<TrackClusterMatchCollection>(tag);
      Print(ih);
    }
  }
}

void 
mu2e::TrackClusterMatchPrinter::Print(const art::Handle<TrackClusterMatchCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::TrackClusterMatchPrinter::Print(const art::ValidHandle<TrackClusterMatchCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::TrackClusterMatchPrinter::Print(const TrackClusterMatchCollection& coll, std::ostream& os) {
  if(verbose()<1) return;
  os << "TrackClusterMatchCollection has " << coll.size() << " matches\n";
  if(verbose()==1) PrintListHeader();
  int i = 0;
  for(const auto& obj: coll) Print(obj, i++);
}

void 
mu2e::TrackClusterMatchPrinter::Print(const art::Ptr<TrackClusterMatch>& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;
  Print(*obj,ind);
}

void 
mu2e::TrackClusterMatchPrinter::Print(const mu2e::TrackClusterMatch& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;

  os << std::setiosflags(std::ios::fixed | std::ios::right);
  if(ind>=0) os << std::setw(4) << ind;

  if(verbose()==1) {
    os       
      << " " << std::setw(8) << std::setprecision(1) << obj.du()
      << " " << std::setw(8) << std::setprecision(1) << obj.dv()
      << " " << std::setw(8) << std::setprecision(3) << obj.ep()
      << " " << std::setw(8) << std::setprecision(3) << obj.chi2()
      << " " << std::setw(8) << std::setprecision(3) << obj.chi2_time()
      << std::endl;
  } else if(verbose()==2) {
    os
      << " icl: " << std::setw(4) << obj.icl()
      << " iex: " << std::setw(4) << obj.iex()
      << "  trk: "
      << " " << std::setw(8) << std::setprecision(1) << obj.xtrk()
      << " " << std::setw(8) << std::setprecision(1) << obj.ytrk()
      << " " << std::setw(8) << std::setprecision(1) << obj.ztrk()
      << " " << std::setw(8) << std::setprecision(1) << obj.ttrk() << "\n";
    os
      << "  nvec: " 
      << " " << std::setw(8) << std::setprecision(1) << obj.nx()
      << " " << std::setw(8) << std::setprecision(1) << obj.ny()
      << " " << std::setw(8) << std::setprecision(1) << obj.nz()
      << "  dx,y,z: " 
      << " " << std::setw(8) << std::setprecision(1) << obj.dx()
      << " " << std::setw(8) << std::setprecision(1) << obj.dy()
      << " " << std::setw(8) << std::setprecision(1) << obj.dz() << "\n";
    os
      << "  du,v: " 
      << " " << std::setw(8) << std::setprecision(1) << obj.du()
      << " " << std::setw(8) << std::setprecision(1) << obj.dv()
      << "  dt: " 
      << " " << std::setw(8) << std::setprecision(1) << obj.dt() << "\n";
    os
      << "  e/p: " << std::setw(8) << std::setprecision(3) << obj.ep()
      << "  chi2: " << std::setw(8) << std::setprecision(3) << obj.chi2()
      << "  chi2_time: " << std::setw(8) << std::setprecision(3) 
         << obj.chi2_time() << "\n";
    os
      << "  int_depth: " << std::setw(8) << std::setprecision(1) << obj.int_depth()
      << "  ds: " << std::setw(8) << std::setprecision(1) << obj.ds()
      << "  dr: " << std::setw(8) << std::setprecision(1) << obj.dr()
      << "  sint: " << std::setw(8) << std::setprecision(4) << obj.sint()
      << std::endl;
    //printf("debug %f drdebug\n",obj.dr());
    //printf("debug %f sidebug\n",obj.sint());
  }

}

void 
mu2e::TrackClusterMatchPrinter::PrintHeader(const std::string& tag, std::ostream& os) {
  if(verbose()<1) return;
  os << "\nProductPrint " << tag << "\n";
}

void 
mu2e::TrackClusterMatchPrinter::PrintListHeader(std::ostream& os) {
  if(verbose()<1) return;
  os << "ind       du       dv      e/p      chi2   chi2_time\n";

}

