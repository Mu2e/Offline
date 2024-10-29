
#include "Offline/Print/inc/StepPointMCPrinter.hh"
#include "art/Framework/Principal/Provenance.h"
#include <iomanip>
#include <string>

void mu2e::StepPointMCPrinter::Print(art::Event const& event,
                                     std::ostream& os) {
  if (verbose() < 1) return;
  if (tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector<art::Handle<StepPointMCCollection> > vah =
        event.getMany<StepPointMCCollection>();
    for (auto const& ah : vah) Print(ah);
  } else {
    // print requested instances
    for (const auto& tag : tags()) {
      auto ih = event.getValidHandle<StepPointMCCollection>(tag);
      Print(ih);
    }
  }
}

void mu2e::StepPointMCPrinter::Print(
    const art::Handle<StepPointMCCollection>& handle, std::ostream& os) {
  if (verbose() < 1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back();  // remove trailing dot
  PrintHeader(tag, os);
  Print(*handle);
}

void mu2e::StepPointMCPrinter::Print(
    const art::ValidHandle<StepPointMCCollection>& handle, std::ostream& os) {
  if (verbose() < 1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back();  // remove trailing dot
  PrintHeader(tag, os);
  Print(*handle);
}

void mu2e::StepPointMCPrinter::Print(const StepPointMCCollection& coll,
                                     std::ostream& os) {
  if (verbose() < 1) return;
  os << "StepPointMCCollection has " << coll.size() << " steps\n";
  if (verbose() == 1) PrintListHeader();
  int i = 0;
  for (const auto& obj : coll) Print(obj, i++);
}

void mu2e::StepPointMCPrinter::Print(const art::Ptr<StepPointMC>& obj, int ind,
                                     std::ostream& os) {
  if (verbose() < 1) return;
  Print(*obj, ind);
}

void mu2e::StepPointMCPrinter::Print(const mu2e::StepPointMC& obj, int ind,
                                     std::ostream& os) {
  if (verbose() < 1) return;

  if (obj.momentum().mag() < _pCut) return;

  int pkey = -1;
  if (obj.simParticle()) pkey = int(obj.simParticle().key());

  os << std::setiosflags(std::ios::fixed | std::ios::right);
  if (ind >= 0) os << std::setw(4) << ind;

  if (verbose() == 1) {
    os << " " << std::setw(8) << pkey << " " << std::setw(7) << obj.volumeId()
       << " " << std::setw(10) << std::setprecision(5) << obj.totalEDep() << " "
       << std::setw(10) << std::setprecision(5) << obj.nonIonizingEDep()
       << "    "
       << " " << std::setw(10) << std::setprecision(3) << obj.position().x()
       << " " << std::setw(10) << std::setprecision(3) << obj.position().y()
       << " " << std::setw(10) << std::setprecision(3) << obj.position().z()
       << " "
       << " " << std::setw(8) << obj.momentum().mag()
       << " " << std::setw(8) << obj.time() << " "
       << " " << std::setw(8) << obj.properTime() << " "
       << " " << std::setiosflags(std::ios::left) << obj.endProcessCode().name()
       << std::endl;

  } else {
    os << "  parentKey: " << std::setw(8) << pkey << "  vol: " << std::setw(5)
       << obj.volumeId() << "  eDep: " << std::setw(8) << std::setprecision(5)
       << obj.totalEDep() << "  nonIonEDep: " << std::setw(8)
       << std::setprecision(5) << obj.nonIonizingEDep() << "\n";
    os << "  pos: " << std::setw(8) << std::setprecision(1)
       << obj.position().x() << " " << std::setw(8) << std::setprecision(1)
       << obj.position().y() << " " << std::setw(8) << std::setprecision(1)
       << obj.position().z() << "  mom: " << std::setw(8)
       << std::setprecision(1) << obj.momentum().x() << " " << std::setw(8)
       << std::setprecision(1) << obj.momentum().y() << " " << std::setw(8)
       << std::setprecision(1) << obj.momentum().z() << "\n";

    os << "  time: " << std::setw(8) << std::setprecision(1) << obj.time()
       << "  ptime: " << std::setw(8) << std::setprecision(1)
       << obj.properTime() << "  step: " << std::setw(8) << std::setprecision(1)
       << obj.stepLength()
       << "    endProc: " << std::setiosflags(std::ios::left)
       << obj.endProcessCode().name() << std::endl;

  }  // end if verbose
}

void mu2e::StepPointMCPrinter::PrintHeader(const std::string& tag,
                                           std::ostream& os) {
  if (verbose() < 1) return;
  os << "\nProductPrint " << tag << std::endl;
}

void mu2e::StepPointMCPrinter::PrintListHeader(std::ostream& os) {
  if (verbose() < 1) return;
  os << "ind     parent     vol     eDep     noIonEDep         Position        "
        "               P      time    properTime endProc"
     << std::endl;
}
