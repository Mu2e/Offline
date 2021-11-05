// vim: set sw=2 expandtab :
#include "art/Framework/Art/mu2eapp.h"
#include "art/Framework/Art/detail/info_success.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include <iostream>
#include <fstream>
#include <cstdlib>

using namespace std;


void mu2eBanner() {
  char* str;
  string art,root,kinkal,rdir,dir,stub,btime,bdir;
  str = getenv("ART_VERSION");
  if(str) art = str;
  str = getenv("KINKAL_VERSION");
  if(str) kinkal = str;
  str = getenv("ROOT_VERSION");
  if(str) root = str;
  str = getenv("MUSE_WORK_DIR");
  if(str) rdir = str;
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

  cout << "   ************************** Mu2e Offline **************************" << endl;
  cout << "     art "<< art <<"    root " << root 
       << "    KinKal " << kinkal << endl;
  cout << "     build  " << dir << endl;
  cout << "     build  " << stub << "    " << btime << endl;
  cout << "   ******************************************************************" << endl;

}



int
main(int argc, char* argv[])
{

  mu2eBanner();

  int result = mu2eapp(argc, argv);
  mf::EndMessageFacility();
  if (result == art::detail::info_success()) {
    return 0;
  }
  cout << "Art has completed and will exit with status " << result << "." << endl;
  return result;
}

// Local Variables:
// mode: c++
// End:
