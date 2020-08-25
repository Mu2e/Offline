
#include "Print/inc/SimStageEfficiencyPrinter.hh"
#include "art/Framework/Principal/Provenance.h"
#include <string>
#include <iomanip>

void
mu2e::SimStageEfficiencyPrinter::PrintRun(art::Run const& run,
                                std::ostream& os) {
  if(verbose()<1) return;
  if(tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector< art::Handle<SimStageEfficiency> > vah;
    run.getManyByType(vah);
    for (auto const & ah : vah) Print(ah, os);
  } else {
    // print requested instances
    for(const auto& tag : tags() ) {
      auto ih = run.getValidHandle<SimStageEfficiency>(tag);
      Print(ih, os);
    }
  }
}

void
mu2e::SimStageEfficiencyPrinter::Print(
             const art::Handle<SimStageEfficiency>& handle,
             std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  if (_forDb) {
    tag = handle.provenance()->moduleLabel();
  }
  PrintHeader(tag,os);
  Print(*handle, os);
}

void
mu2e::SimStageEfficiencyPrinter::Print(
             const art::ValidHandle<SimStageEfficiency>& handle,
             std::ostream& os) {
  if(verbose()<1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back(); // remove trailing dot
  if (_forDb) {
    tag = handle.provenance()->moduleLabel();
  }
  PrintHeader(tag,os);
  Print(*handle,os);
}

void
mu2e::SimStageEfficiencyPrinter::Print(
  const SimStageEfficiency& obj, std::ostream& os) {
  if(verbose()<1) return;
  if (!_forDb) {
    os << "SimStageEfficiency\n";
    os << "Passed Events = " << obj.numerator() << ", All Events = " << obj.denominator() << ", Efficiency = " << obj.efficiency() << std::endl;
  }
  else {
    os << obj.numerator() << "," << obj.denominator() << "," << obj.efficiency() << std::endl;
  }
  //  int i = 0;
  //  for(const auto& obj: coll) Print(obj, i++, os);
}

void
mu2e::SimStageEfficiencyPrinter::PrintHeader(const std::string& tag, std::ostream& os) {
  if(verbose()<1) return;
  if(verbose()==1) {
    if (!_forDb) {
      os << "\nProductPrint " << tag << "\n";
    }
    else {
      os << tag << ",";
    }
  }
}

