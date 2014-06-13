#ifndef TStnEvent_hh
#define TStnEvent_hh

#include "AbsEvent.hh"

class TStnNode;
class TStnDataBlock;
class TStnErrorLogger;

#include "TObjArray.h"

class TStnEvent : public TObject {
public:
  TObjArray*        fListOfInputNodes;	// list of input nodes
  TObjArray*        fListOfOutputNodes;	// list of output nodes
  TObjArray*        fListOfNodes;	// list of registered branches
  int               fLastNumber;	// # of the last printed event
  int               fLastRunNumber;	// run # of the last printed event
  TObjArray*        fListOfObjects;     // 
  TObjArray*        fDropList;          // !
  TStnErrorLogger*  fErrorLogger;	// !
  TObjArray*        fListOfHptl;	// ! list of high-Pt leptons
  Int_t             fCurrentEntry;	// ! current entry in the chain
  Int_t             fCurrentTreeEntry;  // ! current entry in the file
  TObjArray*        fListOfUnusedNodes; // ! list of unused nodes
  Int_t             fEventNumber ;      // !
  Int_t             fRunNumber   ;      // !
  Int_t             fSectionNumber;     // !
//------------------------------------------------------------------------------
//  function members
//------------------------------------------------------------------------------
  TStnEvent();
  virtual ~TStnEvent();

  int    Init(AbsEvent*       event, int mode);

				// for now - from GENP
  int    InitFromGenp();
//-----------------------------------------------------------------------------
//  accessors,  nodes are persistent objects
//-----------------------------------------------------------------------------
  Int_t         EventNumber         () { return fEventNumber; }
  Int_t         RunNumber           () { return fRunNumber; }
  Int_t         SectionNumber       () { return fSectionNumber; }
  TObjArray*    GetListOfNodes      () { return fListOfNodes; }
  TObjArray*    GetListOfInputNodes () { return fListOfInputNodes; }
  TObjArray*    GetListOfOutputNodes() { return fListOfOutputNodes; }
  TObjArray*    GetDropList         () { return fDropList; }
  Int_t         GetCurrentEntry     () { return fCurrentEntry; }
  Int_t         GetCurrentTreeEntry () { return fCurrentTreeEntry; }

  TObjArray*    GetListOfObjects    () { return fListOfObjects; }
  TObjArray*    GetListOfHptl       () { return fListOfHptl; }
  TObjArray*    GetListOfUnusedNodes() { return fListOfUnusedNodes; }

					// read Entry from current tree in the
					// chain

  Int_t         ReadTreeEntry        (Int_t Entry);

				// to pointer to the data block, 
				// corresponding to the branch with the 
				// given name

  TStnDataBlock*   GetDataBlock       (const char* branch_name);

         // get the block, call GetEntry and return the pointer.  
         // If null, then branch does not exist, or user did
         // not register it in user's TStnModule
  TStnDataBlock*   UnpackDataBlock    (const char* branch_name);

				// pointer to pointer to the data block, 
				// corresponding to the branch with the 
				// given name

  TStnDataBlock**  GetDataBlockAddress(const char* branch_name);

  TStnNode*   FindNode(const char* name) { 
    return (TStnNode*) fListOfNodes->FindObject(name);
  }

  virtual TObject*   FindObject(const char* name) const { 
    return fListOfObjects->FindObject(name);
  }

  virtual TObject*   FindObject(const TObject* obj) const { 
    return fListOfObjects->FindObject(obj);
  }

  TStnErrorLogger*  GetErrorLogger() { return fErrorLogger; }
//-----------------------------------------------------------------------------
// setters/modifiers
//-----------------------------------------------------------------------------
  void SetEventNumber(Int_t RunNumber, Int_t EventNumber, Int_t SectionNumber);
  void SetErrorLogger(TStnErrorLogger* Logger) { fErrorLogger  = Logger; }
  void SetCurrentEntry(Int_t             Entry) { fCurrentEntry = Entry ; }
  void SetCurrentTreeEntry(Int_t         Entry) { fCurrentTreeEntry = Entry ; }

  Int_t AddDataBlock(const char* branch_name, 
		     const char* class_name,
		     TStnNode*&  node);
  
  Int_t AddOutputBlock(const char* BranchName, TStnDataBlock* Block);

  void  AddObject(TObject* obj) { fListOfObjects->Add(obj); }
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  void  Clear(Option_t* opt="");

  void  Print(Option_t* opt="") const ;

  ClassDef(TStnEvent,50)
};
#endif
