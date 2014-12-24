//
// Plugin to test some features of FHICL.
//
// $Id: FHICLTest00_module.cc,v 1.4 2013/10/21 21:01:23 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/10/21 21:01:23 $
//
// Original author Rob Kutschke.
//

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"

// Other includes
#include "messagefacility/MessageLogger/MessageLogger.h"

// C++ includes.
//#include <iostream>
#include <string>

using namespace std;

namespace mu2e {

  class FHICLTest00 : public art::EDAnalyzer {
  public:
    explicit FHICLTest00(fhicl::ParameterSet const& pset);
    virtual ~FHICLTest00() {}

    void analyze( art::Event const& e);

  private:

  };

  FHICLTest00::FHICLTest00(fhicl::ParameterSet const& pset):
    art::EDAnalyzer(pset){

    mf::LogVerbatim("Tracing") << "FHICLTest00:analyze: pset:    "
                               << pset.to_string();

    string s(pset.get<string>("p00"));
    mf::LogVerbatim("Tracing") << "FHICLTest00:analyze: p00 as string: " << s;

    int i(pset.get<int>("p00"));
    mf::LogVerbatim("Tracing") << "FHICLTest00:analyze: p00 as int:    " << i;

    double a(pset.get<double>("p00"));
    mf::LogVerbatim("Tracing") << "FHICLTest00:analyze: p00 as double: " << a;

    bool b(pset.get<bool>("p01t"));
    mf::LogVerbatim("Tracing") << "FHICLTest00:analyze: p01t as bool: " << b;

    b = pset.get<bool>("p01f");
    mf::LogVerbatim("Tracing") << "FHICLTest00:analyze: p01f as bool: " << b;

    s = pset.get<string>("p02t");
    mf::LogVerbatim("Tracing") << "FHICLTest00:analyze: p02t as string: " << s;

    s = pset.get<string>("p03t");
    mf::LogVerbatim("Tracing") << "FHICLTest00:analyze: p03t as string: " << s;

    a = pset.get<double>("p04");
    mf::LogVerbatim("Tracing") << "FHICLTest00:analyze: p04 as double: " << a;

    i = pset.get<int>("p04");
    mf::LogVerbatim("Tracing") << "FHICLTest00:analyze: p04 as int:    " << i;

    b = pset.get<bool>("p03t");
    mf::LogVerbatim("Tracing") << "FHICLTest00:analyze: p03t as bool: " << b;

    b = pset.get<bool>("p02t");
    mf::LogVerbatim("Tracing") << "FHICLTest00:analyze: p02t as bool: " << b;

    b = pset.get<bool>("p02f");
    mf::LogVerbatim("Tracing") << "FHICLTest00:analyze: p02f as bool: " << b;

  }

  void FHICLTest00::analyze(art::Event const& ) {
  }

}

using mu2e::FHICLTest00;
DEFINE_ART_MODULE(FHICLTest00)
