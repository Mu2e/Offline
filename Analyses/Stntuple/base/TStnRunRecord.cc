///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "Stntuple/base/TStnRunRecord.hh"

ClassImp(TStnRunRecord)

//-----------------------------------------------------------------------------
TStnRunRecord::TStnRunRecord(int RunNumber) {
  fRunNumber = RunNumber;
}

//-----------------------------------------------------------------------------
TStnRunRecord::~TStnRunRecord() {
}

//-----------------------------------------------------------------------------
int TStnRunRecord::AddEvent(int EventNumber) {

  int  rc;
  int  found = 0;
  int  nev   = fEvent.size();

  for (int i=0; i<nev; i++) {
    if (fEvent[i] == EventNumber) {
      found = 1;
      break;
    }
  }

  if (!found) {
    fEvent.push_back(EventNumber);
    rc = 0;
  }
  else {
    rc = -1;
  }

  return rc;
}

//-----------------------------------------------------------------------------
void TStnRunRecord::Clear(Option_t* Opt) {
  fRunNumber = -1;
  fEvent.clear();
}

//-----------------------------------------------------------------------------
void TStnRunRecord::Print(Option_t* Opt) const {
  printf(" -- Run Number = %10i\n",fRunNumber);
  int k = 0;

  int nev = fEvent.size();

  for (int i=0; i<nev; i++) {
    printf("%10i",fEvent[i]);
    k = k+1;
    if (k == 10) {
      k = 0;
      printf("\n");
    }
  }

  if (k != 0) {
    printf("\n");
  }

  return;
}

