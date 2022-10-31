#include "Offline/Print/inc/STMWaveformDigiPrinter.hh"
#include "art/Framework/Principal/Provenance.h"
#include <iomanip>
#include <string>

void mu2e::STMWaveformDigiPrinter::Print(art::Event const& event,
                                         std::ostream& os) {
  if (verbose() < 1) return;
  if (tags().empty()) {
    // if a list of instances not specified, print all instances
    std::vector<art::Handle<STMWaveformDigiCollection> > vah =
        event.getMany<STMWaveformDigiCollection>();
    for (auto const& ah : vah) Print(ah);
  } else {
    // print requested instances
    for (const auto& tag : tags()) {
      auto ih = event.getValidHandle<STMWaveformDigiCollection>(tag);
      Print(ih);
    }
  }
}

void mu2e::STMWaveformDigiPrinter::Print(
    const art::Handle<STMWaveformDigiCollection>& handle, std::ostream& os) {
  if (verbose() < 1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back();  // remove trailing dot
  PrintHeader(tag, os);
  Print(*handle);
}

void mu2e::STMWaveformDigiPrinter::Print(
    const art::ValidHandle<STMWaveformDigiCollection>& handle,
    std::ostream& os) {
  if (verbose() < 1) return;
  // the product tags with all four fields, with underscores
  std::string tag = handle.provenance()->productDescription().branchName();
  tag.pop_back();  // remove trailing dot
  PrintHeader(tag, os);
  Print(*handle);
}

void mu2e::STMWaveformDigiPrinter::Print(const STMWaveformDigiCollection& coll,
                                         std::ostream& os) {
  if (verbose() < 1) return;
  os << "STMWaveformDigiCollection has " << coll.size() << " entries\n";
  if (verbose() == 1) PrintListHeader();
  int i = 0;
  for (const auto& obj : coll) Print(obj, i++);
}

void mu2e::STMWaveformDigiPrinter::Print(const art::Ptr<STMWaveformDigi>& obj,
                                         int ind, std::ostream& os) {
  if (verbose() < 1) return;
  Print(*obj, ind);
}

void mu2e::STMWaveformDigiPrinter::Print(const mu2e::STMWaveformDigi& obj,
                                         int ind, std::ostream& os) {
  if (verbose() < 1) return;

  os << std::setiosflags(std::ios::fixed | std::ios::right);
  if (ind >= 0) os << std::setw(4) << ind;

  auto const& adcs = obj.adcs();
  os << "  " << std::setw(6) << obj.trigTimeOffset() << std::setw(12)
     << adcs.size() << "  ";
  size_t n = adcs.size();
  if (verbose() <= 1 && n > 5) {
    n = 5;
  }
  for (size_t i = 0; i < n; i++) {
    os << " " << std::setw(6) << adcs[i];
  }
  if (n < adcs.size()) {
    os << " ...";
  }
  os << std::endl;
}

void mu2e::STMWaveformDigiPrinter::PrintHeader(const std::string& tag,
                                               std::ostream& os) {
  if (verbose() < 1) return;
  os << "\nProductPrint " << tag << "\n";
}

void mu2e::STMWaveformDigiPrinter::PrintListHeader(std::ostream& os) {
  if (verbose() < 1) return;
  os << " ind     time      Nadc       Waveform\n";
}
