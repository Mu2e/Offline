//
// Dump information about all data products in the file, including
// event, run and subrun data products.
//
//
// Original author Rob Kutschke
//

// C++ includes.
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Modules/ProvenanceDumper.h"
#include "fhiclcpp/ParameterSet.h"

using namespace std;

namespace mu2e {

  class DataProductDumpDetail;

  typedef art::ProvenanceDumper<DataProductDumpDetail> DataProductDump;

}

class mu2e::DataProductDumpDetail {
public:

  explicit DataProductDumpDetail(fhicl::ParameterSet const& pset);

  typedef std::vector<art::Provenance> provs_type;

  void analyze(const art::Event& e);
  void preProcessEvent() { eventProvenances_.clear(); }
  void processEventProvenance(art::Provenance const & prov) { eventProvenances_.push_back(prov); }
  void postProcessEvent() { process("Event", eventProvenances_); }
  void preProcessSubRun() { subRunProvenances_.clear(); }
  void processSubRunProvenance(art::Provenance const & prov) { subRunProvenances_.push_back(prov); }
  void postProcessSubRun() { process("SubRun", subRunProvenances_); }
  void preProcessRun() { runProvenances_.clear(); }
  void processRunProvenance(art::Provenance const & prov) { runProvenances_.push_back(prov); }
  void postProcessRun() { process("Run", runProvenances_); }

private:
  void process(std::string const & label, provs_type const & provs) const;
  void print(provs_type const & provs) const;

  provs_type eventProvenances_;
  provs_type subRunProvenances_;
  provs_type runProvenances_;
};

mu2e::DataProductDumpDetail::DataProductDumpDetail(fhicl::ParameterSet const& pset)
  :
  eventProvenances_(),
  subRunProvenances_(),
  runProvenances_()
{
}

void mu2e::DataProductDumpDetail::process(std::string const & label,
                                          provs_type const &provs) const {
  std::cout << "Found "
            << provs.size()
            << " data products in this "
            << label
            << std::endl;

  if (provs.size()) {
    std::cout << "Data products: " << std::endl;
    print(provs);
  }
}

// A utility function used by analyze.
void mu2e::DataProductDumpDetail::print( provs_type const& provs) const {

  // Column headings
  string head0("Friendly Class Name");
  string head1("Module Label");
  string head2("Instance Name");
  string head3("Process Name");
  string head4("Product ID");

  // Compute lengths of column headings.
  string::size_type maxlclass(head0.size());
  string::size_type maxlmodule(head1.size());
  string::size_type maxlinstance(head2.size());
  string::size_type maxlprocess(head3.size());
  string::size_type maxlprodId(head4.size());

  // Loop over all products and compute maximum lengths for each column.
  for ( provs_type::const_iterator i=provs.begin();
        i != provs.end(); ++i ){
    art::Provenance const& prov = *i;

    string::size_type lclass, lmodule, linstance, lprocess;
    lclass    = prov.friendlyClassName().size();
    lmodule   = prov.moduleLabel().size();
    linstance = prov.productInstanceName().size();
    lprocess  = prov.processName().size();

    maxlclass    = ( lclass    > maxlclass )    ? lclass    : maxlclass;
    maxlmodule   = ( lmodule   > maxlmodule )   ? lmodule   : maxlmodule;
    maxlinstance = ( linstance > maxlinstance ) ? linstance : maxlinstance;
    maxlprocess  = ( lprocess  > maxlprocess )  ? lprocess  : maxlprocess;
  }

  // Print table header.
  cout << setw(maxlclass)    << head0 << "  "
       << setw(maxlmodule)   << head1 << "  "
       << setw(maxlinstance) << head2 << "  "
       << setw(maxlprocess)  << head3 << "     "
       << setw(maxlprodId)   << head4
       << endl;

  // Printer body of header.
  for ( provs_type::const_iterator i=provs.begin();
        i != provs.end(); ++i ){
    art::Provenance const& prov = *i;
    cout << setw(maxlclass)    << prov.friendlyClassName()   << "  "
         << setw(maxlmodule)   << prov.moduleLabel()         << "  "
         << setw(maxlinstance) << prov.productInstanceName() << "  "
         << setw(maxlprocess)  << prov.processName()         << "  "
         << setw(maxlprodId)   << prov.productID()
         << endl;
  }

  cout << endl;

} // end print

DEFINE_ART_MODULE(mu2e::DataProductDump);
