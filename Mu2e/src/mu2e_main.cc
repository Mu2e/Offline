// vim: set sw=2 expandtab :

// Initialize ExitCodePrinter as early as possible, before other
// include files, to have its destructor called as late as possible.
#include "art/Framework/Art/detail/ExitCodePrinter.h"
namespace {
  art::detail::ExitCodePrinter p;
}

#include "art/Framework/Art/artapp.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <iostream>
#include <fstream>
#include <cstdlib>

using namespace std;


void mu2eBanner() {
  char* strs;
  char* str;
  string art,root,kinkal,rdir,dir,stub,btime,bdir;

  cout << "   ************************** Mu2e Offline **************************" << endl;

  strs = getenv("OFFLINE_BANNER");
  str = getenv("MUSE_WORK_DIR");
  if(str) rdir = str;
  if(strs) { // spack only build
    cout << "     " << strs << endl;
  } else if(str) { // muse build
    str = getenv("ART_VERSION");
    if(str) art = str;
    str = getenv("KINKAL_VERSION");
    if(str) kinkal = str;
    str = getenv("ROOT_VERSION");
    if(str) root = str;
    str = getenv("MUSE_STUB");
    if(str) stub = str;

    btime = "musebuild file not found";
    dir = rdir;
    ifstream mfile(rdir+"/build/"+stub+"/.musebuild");
    if (mfile.is_open()) {
      getline(mfile,btime);
      if(btime.size()>17) btime = btime.substr(0,17);
      getline(mfile,bdir);
      if(bdir.size()>0) dir = bdir;
    }
    mfile.close();

    cout << "     art "<< art <<"    root " << root
         << "    KinKal " << kinkal << endl;
    cout << "     build  " << dir << endl;
    cout << "     build  " << stub << "    " << btime << endl;
  } // end muse format
  cout << "   ******************************************************************" << endl;

}

int
main(int argc, char* argv[])
{
  mu2eBanner();
  p = artapp(argc, argv, false);
  mf::EndMessageFacility();
  return p.exitcode();
}

// Local Variables:
// mode: c++
// End:
