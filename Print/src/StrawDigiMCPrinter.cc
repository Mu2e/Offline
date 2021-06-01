
#include "Print/inc/StrawDigiMCPrinter.hh"
#include "art/Framework/Principal/Provenance.h"
#include <string>
#include <iomanip>

void 
mu2e::StrawDigiMCPrinter::Print(art::Event const& event,
				std::ostream& os) {
  if(verbose()<1) return;
  if(tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector< art::Handle<StrawDigiMCCollection> > vah = event.getMany<StrawDigiMCCollection>();
    for (auto const & ah : vah) Print(ah);
  } else {
    // print requested instances
    for(const auto& tag : tags() ) {
      auto ih = event.getValidHandle<StrawDigiMCCollection>(tag);
      Print(ih);
    }
  }
}

void 
mu2e::StrawDigiMCPrinter::Print(const 
				art::Handle<StrawDigiMCCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::StrawDigiMCPrinter::Print(const art::ValidHandle<StrawDigiMCCollection>& handle,
				std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  PrintHeader(tag,os);
  Print(*handle);
}

void 
mu2e::StrawDigiMCPrinter::Print(const StrawDigiMCCollection& coll, std::ostream& os) {
  if(verbose()<1) return;
  os << "StrawDigiMCCollection has " << coll.size() << " hits\n";
  if(verbose()==1) PrintListHeader();
  int i = 0;
  for(const auto& obj: coll) Print(obj, i++);
}

void 
mu2e::StrawDigiMCPrinter::Print(const art::Ptr<StrawDigiMC>& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;
  Print(*obj,ind);
}

void 
mu2e::StrawDigiMCPrinter::Print(const mu2e::StrawDigiMC& obj, int ind, std::ostream& os) {
  if(verbose()<1) return;

  os << std::setiosflags(std::ios::fixed | std::ios::right);
  if(ind>=0) os << std::setw(4) << ind;

  // energy accessor is not protected against bad Ptr
  double energy = 0.0, tenergy0 = 0.0, tenergy1 = 0.0;
  bool eOK = true;
  if(!obj.strawGasStep(StrawEnd::cal).isAvailable() ||
      !obj.strawGasStep(StrawEnd::hv).isAvailable() ||
      !obj.strawGasStep(StrawEnd::cal).isAvailable() ||
      !obj.strawGasStep(StrawEnd::hv).isAvailable() ) eOK = false;
  if(eOK) {
    energy = obj.energySum();
    tenergy0 = obj.triggerEnergySum(StrawEnd::cal);
    tenergy1 = obj.triggerEnergySum(StrawEnd::hv);
  }

  // now check if the basic StepPointMC's are also available
  auto const& a0 = obj.strawGasStep(StrawEnd::cal);
  auto const& a1 = obj.strawGasStep(StrawEnd::hv);

  if(!(a0.isAvailable() && a1.isAvailable()) ) {
    // give up at this point since basically no accessors work
    std::cout << " StepPointMC products not available" << std::endl;
    return;
  }

  if (verbose()==1) {
    os << " " << std::setw(5) << obj.strawId().asUint16()
       << " " 
       << " " << std::setw(8) << std::setprecision(2) 
           << obj.wireEndTime(StrawEnd::cal)
       << " " << std::setw(8) << std::setprecision(2) 
           << obj.wireEndTime(StrawEnd::hv)
       << " " << std::setw(8) << std::setprecision(4) << energy
       << " " << std::setw(6) << obj.isCrossTalk(StrawEnd::cal)
       << " " << std::setw(6) << obj.isCrossTalk(StrawEnd::hv)
       << " " ;
    os << std::endl;
  } else if (verbose()==2) {

    os << " index :" << std::setw(5) << obj.strawId().asUint16()
       << " " 
       << " time0: " << std::setw(8) << std::setprecision(2) 
                  << obj.wireEndTime(StrawEnd::cal)
       << " time1: " << std::setw(8) << std::setprecision(2) 
                  << obj.wireEndTime(StrawEnd::hv)
       << " " ;
    os << std::endl;
    os << " trigEnergy0: " << std::setw(8) << std::setprecision(4) 
              << tenergy0
       << " trigEnergy1: " << std::setw(8) << std::setprecision(4) 
              << tenergy1
       << " cross0: " << std::setw(2) << obj.isCrossTalk(StrawEnd::cal)
       << " cross1: " << std::setw(2) << obj.isCrossTalk(StrawEnd::hv)
       << " " ;
    os << std::endl;
    CLHEP::HepLorentzVector const& v0 = obj.clusterPosition(StrawEnd::cal);
    os << " position0: " 
       << std::setw(10) << std::setprecision(4) << v0.x()
       << std::setw(10) << std::setprecision(4) << v0.y()
       << std::setw(10) << std::setprecision(4) << v0.z()
       << std::setw(10) << std::setprecision(4) << v0.t()
       << " " ;
    os << std::endl;
    CLHEP::HepLorentzVector const& v1 = obj.clusterPosition(StrawEnd::hv);
    os << " position1: " 
       << std::setw(10) << std::setprecision(4) << v1.x()
       << std::setw(10) << std::setprecision(4) << v1.y()
       << std::setw(10) << std::setprecision(4) << v1.z()
       << std::setw(10) << std::setprecision(4) << v1.t()
       << " " ;
    os << std::endl;

    if(!a0.isAvailable()) {
      os << " StepPointMC_0: is invalid " << std::endl;
    } else {
      int isim = -1;
      if(a0->simParticle()) isim = int(a0->simParticle().key());
      os << " StepPointMC_0: " 
	 << "   key: " << std::setw(5) << a0.key()
	 << "   energy: " << std::setw(8) << std::setprecision(6)  
	       << a0->eDep()
	 << "   Simparticle: " << std::setw(5) << isim;
      os << std::endl;
    }

    if(!a1.isAvailable()) {
      os << " StepPointMC_1: is invalid " << std::endl;
    } else {
      int isim = -1;
      if(a1->simParticle()) isim = int(a1->simParticle().key());
      os << " StepPointMC_1: " 
	 << "   key: " << std::setw(5) << a1.key()
	 << "   energy: " << std::setw(8) << std::setprecision(4)  
	       << a1->eDep()
	 << "   SimParticle: " << std::setw(5) << isim;
      os << std::endl;
    }
  } // verbose=2

}

void 
mu2e::StrawDigiMCPrinter::PrintHeader(const std::string& tag, std::ostream& os) {
  if(verbose()<1) return;
  os << "\nProductPrint " << tag << "\n";
}

void 
mu2e::StrawDigiMCPrinter::PrintListHeader(std::ostream& os) {
  if(verbose()<1) return;
  os << " ind StrwInd   time0   time1    drift0   drift1   energy    crossTalk\n";

}

