//-----------------------------------------------------------------------------
//  Dec 15 1999 P.Murat: 
//  --------------------
//  an example showing how to copy a subset of events from the input ntuple 
//  (UC STNTUPLE 5.x converted to ROOT) to output ntuple in the same format
//
//  the user is supposed to supply a function event_is_ok(THBookStnEvent*)
//  returning 1 if the input event has to be copied into output ntuple and
//  0 otherwise
//-----------------------------------------------------------------------------

#ifndef __CINT__
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TInterpreter.h"
#include <Stntuple/THBookStnEvent.hh>
#include <Stntuple/THBookStntuple.hh>
#include <Stntuple/TStntupleEvent.hh>
#include <Stntuple/TStntuple.hh>
#endif

THBookStnEvent* input_event;

//-----------------------------------------------------------------------------
// return 1 if there is at least 1 loose tau candidate
//-----------------------------------------------------------------------------
void event_is_ok () {
  int rc = 0;
  StnTausBlock_t* blk = input_event->fTausBlock;
  for (int i=0; i<blk->Ntau; i++) {
    if (blk->Tstat[i] > 2000.) {
      rc = 1;
      break;
    }
  }
  input_event->SetPassed(rc);
}

//-----------------------------------------------------------------------------
//  `dir' - name of the directory with root files, I'm assuming that all the
// *.root files in `dir' produce one chain
//-----------------------------------------------------------------------------
int copy_hb(const char*  dir_name, 
	    const char*  output_file_name, 
	    int          first_event=0, 
	    int          last_event=-1) 
{
  input_event  = new THBookStnEvent();
				        // define input chain

  THBookStntuple* input = new THBookStntuple("/STNTUPLE/h1",input_event);
  input->Init(dir_name);

  input->CopyChain(output_file_name,"event_is_ok");
  return 0;
}

