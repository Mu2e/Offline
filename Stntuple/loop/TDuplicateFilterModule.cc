///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "Stntuple/base/TStnRunRecord.hh"
#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/obj/TStnHeaderBlock.hh"

#include "Stntuple/loop/TDuplicateFilterModule.hh"

ClassImp(TDuplicateFilterModule)
//_____________________________________________________________________________
TDuplicateFilterModule::TDuplicateFilterModule(const char* name, const char* title):
  TStnModule(name,title) 
{
//-----------------------------------------------------------------------------
// by default print run-level catalog of each file (fPrintLevel=2)
//-----------------------------------------------------------------------------
  fListOfRunRecords = new TObjArray();
  fNDuplicateEvents = 0;
  fCurrentRunRecord = 0;
}

//_____________________________________________________________________________
TDuplicateFilterModule::~TDuplicateFilterModule() {
  delete fListOfRunRecords;
}

//_____________________________________________________________________________
int TDuplicateFilterModule::BeginJob() {
  TStnModule::BeginJob();
  fNDuplicateEvents = 0;
  fCurrentRunRecord = 0;
  return 0;
}


//_____________________________________________________________________________
int TDuplicateFilterModule::BeginRun() {
//---------------------------------------------------------------------------
// Figure out whether we know this run, otherwise install a new record
// Try to find this run in our filing system
//---------------------------------------------------------------------------
  TStnRunRecord* rr;

  int  rn    = GetHeaderBlock()->RunNumber();
  bool found = false;

  int nruns = fListOfRunRecords->GetEntries();

  for (int i=0; i<nruns; i++) {
    rr = (TStnRunRecord*) fListOfRunRecords->At(i);
    if (rr->RunNumber() == rn) {
      fCurrentRunRecord = rr;
      found = true;
    }
  }
				// This run was not yet found, create new one
  if (! found) {
    fCurrentRunRecord = new TStnRunRecord(rn);
    fListOfRunRecords->Add(fCurrentRunRecord);
  }

  return 0;
}

//_____________________________________________________________________________
int TDuplicateFilterModule::Event(int ientry) {

  int passed(1);

  TStnHeaderBlock* header = GetHeaderBlock();

  //  int rn = header->RunNumber    ();
  int ev = header->EventNumber  ();
  //  int rs = header->SectionNumber();

				// current run record is already set

  int rc = fCurrentRunRecord->AddEvent(ev);

  if (rc == -1) {
				// this event is a duplicate
    passed = 0;
    fNDuplicateEvents += 1;

    header->Print("TDuplicateFilterModule::Event WARNING: Duplicate Event");
  }

  SetPassed(passed);

  return 0;
}

//_____________________________________________________________________________
int TDuplicateFilterModule::EndJob() {
  printf(" TDuplicateFilterModule: N(duplicate events) = %10i\n",
	 fNDuplicateEvents);
  return 0;
}
