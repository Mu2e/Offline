#ifndef TStnRunRecord_hh
#define TStnRunRecord_hh
///////////////////////////////////////////////////////////////////////////////
// 2010-01-21 P.Murat: just book-kkeping, to get rid of event duplicates
///////////////////////////////////////////////////////////////////////////////
#include <vector>
#include "TObject.h"
#include "TObjArray.h"

using std::vector;

class TStnRunRecord : public TObject {
public:
  int           fRunNumber;
  vector<int>   fEvent;

   TStnRunRecord(int RunNumber = -1);
  ~TStnRunRecord();

  void   Print(Option_t* Option = "") const ;

  int    RunNumber() { return fRunNumber; }
  int    NEvents  () { return fEvent.size(); }

				// returns  0 if a new event added
				//         -1 if a duplicate

  int    AddEvent (int EventNumber);

  void   Clear(Option_t* Option = "");

  ClassDef(TStnRunRecord,0)
};

#endif
