//
// A plugin to test FileInPath.
//
// $Id: FileInPathTest_module.cc,v 1.4 2011/05/17 22:50:55 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2011/05/17 22:50:55 $
//
// Original author Rob Kutschke
//

// C++ includes.
#include <iostream>
#include <string>

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h"
//#include "art/ParameterSet/FileInPath.h"
#include "art/Framework/Services/Optional/TFileService.h"

using namespace std;

namespace mu2e {

  class Straw;

  class FileInPathTest : public art::EDAnalyzer {
  public:
    
    explicit FileInPathTest(fhicl::ParameterSet const& pset);
    virtual ~FileInPathTest() { }

    void analyze(const art::Event& e);

  private:

  };

  FileInPathTest::FileInPathTest(fhicl::ParameterSet const& pset){
    cout  << "FipTest::constructor." << endl;
    cout  << "This is no longer needed in the art context." << endl;

    // TODO: If FileInPath goes away, then this goes away provided
    //       proper docs exist elsewhere for its replacement.

    //art::FileInPath fp = pset.get<art::FileInPath>("inputfile");
    //art::FileInPath fp2("Mu2eG4/test/ttracker_v0.txt");
    //art::FileInPath fp3("BFieldMaps/MECO/dsmap_unfmt_rad100.dat");
    //cout << "fp is:  " << fp.location()  << " " << fp.isLocal()  << " " << fp.fullPath()  << endl;
    //cout << "fp2 is: " << fp2.location() << " " << fp2.isLocal() << " " << fp2.fullPath() << endl;
    //cout << "fp3 is: " << fp3.location() << " " << fp3.isLocal() << " " << fp3.fullPath() << endl;


  }
  
  void FileInPathTest::analyze(const art::Event& event) {
  }

}  // end namespace mu2e

using mu2e::FileInPathTest;
DEFINE_ART_MODULE(FileInPathTest);
