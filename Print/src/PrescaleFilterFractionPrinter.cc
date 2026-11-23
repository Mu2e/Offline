#include "Offline/Print/inc/PrescaleFilterFractionPrinter.hh"
#include "art/Framework/Principal/Provenance.h"
#include "art/Framework/Principal/SubRun.h"
#include <boost/io/ios_state.hpp>
#include <iomanip>
#include <string>
#include <vector>

void mu2e::PrescaleFilterFractionPrinter::Print(art::Event const& event, std::ostream& os) {}

void mu2e::PrescaleFilterFractionPrinter::PrintSubRun(art::SubRun const& subrun, std::ostream& os) {
  if (verbose() < 1) return;
  if (tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector<art::Handle<PrescaleFilterFraction> > ffl = subrun.getMany<PrescaleFilterFraction>();
    for (auto const& ff : ffl) Print(ff);
  } else {
    // print requested instances
    for (const auto& tag : tags()) {
      auto ff = subrun.getValidHandle<PrescaleFilterFraction>(tag);
      Print(ff);
    }
  }
}

void mu2e::PrescaleFilterFractionPrinter::Print(
    const art::Handle<PrescaleFilterFraction>& handle, std::ostream& os) {
  if (verbose() < 1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back();  // remove trailing dot
  PrintHeader(tag, os);
  Print(*handle);
}

void mu2e::PrescaleFilterFractionPrinter::Print(
    const art::ValidHandle<PrescaleFilterFraction>& handle, std::ostream& os) {
  if (verbose() < 1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back();  // remove trailing dot
  PrintHeader(tag, os);
  Print(*handle);
}

void mu2e::PrescaleFilterFractionPrinter::Print(const mu2e::PrescaleFilterFraction& obj,
    int ind, std::ostream& os) {
  boost::io::ios_all_saver saver(os);  // leave the caller's stream formatting unchanged
  // The first part of the line keeps its original text, which downstream parsers match.
  // The prescale fraction uses significant digits so that small fractions (e.g. 4e-05) stay visible.
  os << std::setiosflags(std::ios::right);
  os << " Fraction passing filter " << std::fixed << std::setprecision(6) << obj.filterFraction()
     << " N Seen " << obj.nSeen()
     << " with prescale fraction " << std::defaultfloat << std::setprecision(6) << obj.prescaleFraction();
  os << " N Passed " << obj.nPassed();
  // a fractional prescale has no integer factor; a discrete one passes 1 in prescale() events
  if (obj.isFractional())
    os << " fractional";
  else if (obj.prescale() == 0)
    os << " discrete disabled";
  else
    os << " discrete 1/" << obj.prescale();
  if(obj.chained())os << " chained";
  os << std::endl;
}

void mu2e::PrescaleFilterFractionPrinter::PrintHeader(const std::string& tag,
    std::ostream& os) {
  if (verbose() < 1) return;
  os << "\nProductPrint " << tag << "\n";
}
