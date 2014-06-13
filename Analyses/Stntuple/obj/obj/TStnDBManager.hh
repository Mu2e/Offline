#ifndef STNTUPLE_TStnDBManager_hh
#define STNTUPLE_TStnDBManager_hh
//-----------------------------------------------------------------------------
// Implementation of storage that is not event related but rather used
// on a run by run basis. It resembles a database with its only
// accessor:
//
//   GetTable("TableName")
//
// In reality the implementation is just simply TObjects sitting in
// the ntuple which are accessed by name.
//
// Example: 
//
//   TStnTriggerTable* t =
//     (TStnTriggerTable*) TStnDBManager::Instance()->GetTable("TriggerTable");
//   if (t == 0)
//      cout << "data not found.\n";
//-----------------------------------------------------------------------------

#include "TObject.h"
#include "TList.h"

class TStnDBManager: public TObject {
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
protected:
  static TStnDBManager*  fgInstance;

  // TList of pointers to Beamline, triggertable etc
  TList * fListOfDbObjects;

  class  Cleaner {
  public: 
    Cleaner ();
    ~Cleaner();
  };
  friend class Cleaner;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:

  TStnDBManager();
  virtual ~TStnDBManager();

  static TStnDBManager*  Instance();
					// ****** init methods

					// ****** accessors
  TObject*  GetTable  (const char* Name);

					// ****** overloaded methods of TObject
  Int_t Read (const char* Name);
  Int_t Write(const char* Name=0, Int_t option=0, Int_t bufsize=0);
  //  void  MakeAvgBeamline();

  ClassDef(TStnDBManager,2)
};
#endif
