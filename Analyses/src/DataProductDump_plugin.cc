//
// Dump information about all data products in the file.
//
// $Id: DataProductDump_plugin.cc,v 1.3 2011/05/17 15:46:09 greenc Exp $
// $Author: greenc $ 
// $Date: 2011/05/17 15:46:09 $
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
#include "art/Framework/Core/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h"

using namespace std;

namespace mu2e {

  class DataProductDump : public art::EDAnalyzer {
  public:
    
    explicit DataProductDump(fhicl::ParameterSet const& pset);
    virtual ~DataProductDump() { }

    typedef std::vector<art::Provenance const*> provs_type;

    void analyze(const art::Event& e, art::EventSetup const&);

    // A utility function to print info about a vector of provenances.
    void print( provs_type const& provs);

  private:

  };

  DataProductDump::DataProductDump(fhicl::ParameterSet const& pset){
  }

  void DataProductDump::analyze(const art::Event& event, art::EventSetup const&) {

    art::Run const& run(event.getRun());
    art::SubRun const& lumi(event.getSubRun());

    // Get provenance information for all data products.
    provs_type evtProvs;
    provs_type runProvs;
    provs_type lumiProvs;
    event.getAllProvenance(evtProvs);
    run.getAllProvenance(runProvs);
    lumi.getAllProvenance(lumiProvs);

    // Print a header for this event.
    cout << "\n\nData producs for: " << event.id() << endl;
    cout << "Data products in this event:      " << evtProvs.size()  << endl;
    cout << "Data products in this run:        " << runProvs.size()  << endl;
    cout << "Data products in this subRun:     " << lumiProvs.size() << endl;
    cout << endl;

    // Print information about each set of data products.
    if ( evtProvs.size() > 0 ){
      cout << "Event data products: " << endl;
      print (evtProvs);
    }

    if ( runProvs.size() > 0 ) {
      cout << "Run data products: " << endl;
      print (runProvs);
    }

    if ( lumiProvs.size() > 0 ) {
      cout << "SubRun data products: " << endl;
      print (lumiProvs);
    }

  } // end analyze

  // A utility function used by analyze.
  void DataProductDump::print( provs_type const& provs){

    // Column headings
    string head0("Friendly Class Name");
    string head1("Module Label");
    string head2("Instance Name");
    string head3("Process Name");

    // Compute lengths of column headings.
    string::size_type maxlclass(head0.size());
    string::size_type maxlmodule(head1.size());
    string::size_type maxlinstance(head2.size());
    string::size_type maxlprocess(head3.size());

    // Loop over all products and compute maximum lengths for each column.
    for ( provs_type::const_iterator i=provs.begin();
          i != provs.end(); ++i ){
      art::Provenance const& prov = **i;

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
         << setw(maxlprocess)  << head3
         << endl;

    // Printer body of header.
    for ( provs_type::const_iterator i=provs.begin();
          i != provs.end(); ++i ){
      art::Provenance const& prov = **i;
      cout << setw(maxlclass)    << prov.friendlyClassName()   << "  "
           << setw(maxlmodule)   << prov.moduleLabel()         << "  "
           << setw(maxlinstance) << prov.productInstanceName() << "  "
           << setw(maxlprocess)  << prov.processName()
           << endl;
    }

    cout << endl;

  } // end print

}  // end namespace mu2e

using mu2e::DataProductDump;
DEFINE_ART_MODULE(DataProductDump);
