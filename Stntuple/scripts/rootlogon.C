//------------------------------------------------------------------------------
//  rootlogon.C: a sample ROOT logon macro allowing use of ROOT script 
//               compiler in CDF RunII environment. The name of this macro file
//               is defined by the .rootrc file
//
//  USESHLIBS variable has to be set to build Stntuple libraries locally:  
//
//  setenv USESHLIBS 1
//
//  Feb 12 2001 P.Murat
//------------------------------------------------------------------------------
{
#include <time.h>
                                // the line below tells CINT where to look for 
				// the include files

  gInterpreter->AddIncludePath(Form("%s/include",
				    gSystem->Getenv("SRT_LOCAL")));

  gInterpreter->AddIncludePath(Form("%s/include",
				    gSystem->Getenv("CDFSOFT2_DIR")));

  gInterpreter->AddIncludePath(Form("%s/tex/cdfnotes",
				    gSystem->Getenv("HOME")));
//   gSystem->SetMakeSharedLib("cd $BuildDir ; g++ -c -g $Opt -pipe -m32 -Wall -W -Woverloaded-virtual -fPIC -pthread $IncludePath $SourceFiles ;  g++ -g $ObjectFiles -shared -Wl,-soname,$LibName.so -m32 $LinkedLibs -o $SharedLib");
//-----------------------------------------------------------------------------
// load in ROOT physics vectors and event generator libraries
//-----------------------------------------------------------------------------
  gSystem->Load("$ROOTSYS/lib/libPhysics.so");
  gSystem->Load("$ROOTSYS/lib/libEG.so");
  gSystem->Load("$ROOTSYS/lib/libMinuit.so");
  gSystem->Load("$ROOTSYS/lib/libFumili.so");
  gSystem->Load("$ROOTSYS/lib/libTree.so");
//-----------------------------------------------------------------------------
//  check batch mode
//-----------------------------------------------------------------------------
  int batch_mode = 0;
  const char* opt;
  int nargs = gApplication->Argc();
  for (int i=1; i<nargs; i++) {
    opt = gApplication->Argv(i);
    if (strcmp(opt,"-b") == 0) {
      batch_mode = 1;
      break;
    }
  }

  printf("   batch_mode = %i\n",batch_mode);
//-----------------------------------------------------------------------------
// STNTUPLE shared libraries are assumed to be built in the private test 
// release area with USESHLIBS environment variable set 
// we always need libStntuple_loop, but the other 2 libs should be loaded in 
// only if we're running bare root
//-----------------------------------------------------------------------------
  const char* exec_name = gApplication->Argv(0);

  if (strstr(exec_name,"root.exe") != 0) {
//-----------------------------------------------------------------------------
// assume STNTUPLE  analysis job
//-----------------------------------------------------------------------------
    if (batch_mode == 1) gSystem->Load("$ROOTSYS/lib/libGui.so");

    gSystem->Load("$MU2E_TEST_RELEASE/lib/libStntuple_base.so");
    gSystem->Load("$MU2E_TEST_RELEASE/lib/libStntuple_obj.so");
    gSystem->Load("$MU2E_TEST_RELEASE/lib/libStntuple_loop.so");
    gSystem->Load("$MU2E_TEST_RELEASE/lib/libStntuple_alg.so");
    gSystem->Load("$MU2E_TEST_RELEASE/lib/libStntuple_val.so");
    gSystem->Load("$MU2E_TEST_RELEASE/lib/libStntuple_ana.so");
  }
                                // print overflows/underflows in the stat box
  gStyle->SetOptStat(11111111);
                                // print fit results in the stat box
  gStyle->SetOptFit(1110);
//-----------------------------------------------------------------------------
//  databases
//-----------------------------------------------------------------------------
//   gSystem->Load("libStntuple_oracle.so");
//-----------------------------------------------------------------------------
//  stnana, if exists and if executable is not mu2e
//-----------------------------------------------------------------------------
//   if (strstr(exec_name,"mu2e.exe") == 0) {
//     if (gSystem->Exec("ls Stntuple/scripts/stnana.C > /dev/null ") == 0) {
//       gInterpreter->LoadMacro("Stntuple/scripts/stnana.C");
//     }
//   }
//-----------------------------------------------------------------------------
// this line reports the process ID which simplifies debugging
//-----------------------------------------------------------------------------
  printf(" process ID: %i\n",gSystem->GetPid());
  TAuthenticate::SetGlobalUser(gSystem->Getenv("USER"));
  gInterpreter->ProcessLine(".! ps | grep root");
}
