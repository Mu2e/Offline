//
// A plugin to test FileInPath.
//
// $Id: FileInPathTest_plugin.cc,v 1.1 2010/09/27 19:45:21 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/09/27 19:45:21 $
//
// Original author Rob Kutschke
//

// C++ includes.
#include <iostream>
#include <string>

// Framework includes.
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/Services/interface/TFileService.h"

using namespace std;

namespace mu2e {

  class Straw;

  class FileInPathTest : public edm::EDAnalyzer {
  public:
    
    explicit FileInPathTest(edm::ParameterSet const& pset);
    virtual ~FileInPathTest() { }

    void analyze(const edm::Event& e, edm::EventSetup const&);

  private:

  };

  FileInPathTest::FileInPathTest(edm::ParameterSet const& pset){
    cout  << "FipTest::constructor." << endl;
    edm::FileInPath fp = pset.getParameter<edm::FileInPath>("inputfile");
    edm::FileInPath fp2("Mu2eG4/test/ttracker_v0.txt");
    edm::FileInPath fp3("BFieldMaps/MECO/dsmap_unfmt_rad100.dat");
    cout << "fp is:  " << fp.location()  << " " << fp.isLocal()  << " " << fp.fullPath()  << endl;
    cout << "fp2 is: " << fp2.location() << " " << fp2.isLocal() << " " << fp2.fullPath() << endl;
    cout << "fp3 is: " << fp3.location() << " " << fp3.isLocal() << " " << fp3.fullPath() << endl;


  }
  
  void FileInPathTest::analyze(const edm::Event& event, edm::EventSetup const&) {
  }

}  // end namespace mu2e

using mu2e::FileInPathTest;
DEFINE_FWK_MODULE(FileInPathTest);
