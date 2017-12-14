
#include "Print/inc/StereoHitPrinter.hh"
#include "art/Framework/Principal/Provenance.h"
#include <string>
#include <iomanip>

void 
mu2e::StereoHitPrinter::Print(art::Event const& event,
				std::ostream& os) {
  if(verbose()<1) return;
  if(_tags.empty()) {
    // if a list of instances not specified, print all instances
    std::vector< art::Handle<StereoHitCollection> > vah;
    event.getManyByType(vah);
    for (auto const & ah : vah) Print(ah);
  } else {
    // print requested instances
    for(const auto& tag : tags() ) {
      auto ih = event.getValidHandle<StereoHitCollection>(tag);
      Print(ih);
    }
  }
}

void 
mu2e::StereoHitPrinter::Print(const art::Handle<StereoHitCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::StereoHitPrinter::Print(const art::ValidHandle<StereoHitCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::StereoHitPrinter::Print(const StereoHitCollection& coll, std::ostream& os) {
  if(verbose()<1) return;
  os << "StereoHitCollection has " << coll.size() << " hits\n";
  if(verbose()==1) PrintListHeader();
  int i = 0;
  for(const auto& obj: coll) Print(obj, i++);
}

void 
mu2e::StereoHitPrinter::Print(const art::Ptr<StereoHit>& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;
  Print(*obj,ind);
}

void 
mu2e::StereoHitPrinter::Print(const mu2e::StereoHit& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;

  os << std::setiosflags(std::ios::fixed | std::ios::right);
  if(ind>=0) os << std::setw(4) << ind;

  os 
    << " " << std::setw(5) << obj.hitIndex1()
    << " " << std::setw(5) << obj.hitIndex2()
    << " " 
    << " " << std::setw(8) << std::setprecision(3) << obj.pos().x()
    << " " << std::setw(8) << std::setprecision(3) << obj.pos().y()
    << " " << std::setw(9) << std::setprecision(3) << obj.pos().z()
    << " " << std::setw(7) << std::setprecision(1) << obj.time()
    << " " << std::setw(8) << std::setprecision(5) << obj.energy()
    << " " << std::setw(6) << std::setprecision(2) << obj.chisq()
    << " " << std::setw(6) << std::setprecision(3) << obj.mvaout()
    << std::endl;

}

void 
mu2e::StereoHitPrinter::PrintHeader(const std::string& tag, std::ostream& os) {
  if(verbose()<1) return;
  os << "\nProductPrint " << tag << "\n";
}

void 
mu2e::StereoHitPrinter::PrintListHeader(std::ostream& os) {
  if(verbose()<1) return;
  os <<  "ind   str1   str2    x        y         z        t        E     chi2   mva\n";

}

void 
mu2e::StereoHitPrinter::set(const fhicl::ParameterSet& pset) {

  fhicl::ParameterSet localPset = 
    pset.get<fhicl::ParameterSet>("StereoHitPrinter",fhicl::ParameterSet());

  setVerbose( localPset.get<int>("verbose",verbose()) );
  _tags = vecstr( localPset.get<vecstr>("inputTags",vecstr()) );

}

