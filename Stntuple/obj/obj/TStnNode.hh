#ifndef STNTUPLE_TStnNode
#define STNTUPLE_TStnNode
//-----------------------------------------------------------------------------
//  definition of the STNTUPLE node (holder for the data block)
//  Author:    Pasha Murat (CDF/FNAL)
//  Date:      Nov 10 2000
// 
//-----------------------------------------------------------------------------

#include "TNamed.h"
class TBranch;
class TStnDataBlock;
class TStnEvent;

class TStnNode: public TNamed {
protected:
//------------------------------------------------------------------------------
//  data members
//------------------------------------------------------------------------------
  TStnDataBlock*  fObject;		// object, corresponding to the branch
  TBranch*        fBranch;
  TStnEvent*      fEvent;
  Int_t           (*fFunc)(TStnDataBlock *, TStnEvent *, Int_t);
  Int_t           fDeleteObject;	// ! 
//------------------------------------------------------------------------------
//  functions
//------------------------------------------------------------------------------
public:
					// ****** constructors and destructor
  TStnNode();
  TStnNode(const char* BranchName, 
	   TClass*     cl,
	   TStnEvent*  Event,
	   Int_t       (*F)(TStnDataBlock *, TStnEvent *, Int_t) = 0);

  TStnNode(TBranch*    Branch, 
	   TClass*     cl,
	   TStnEvent*  Event,
	   Int_t       (*F)(TStnDataBlock *, TStnEvent *, Int_t) = 0);

  TStnNode(const char*    BranchName,
	   TStnDataBlock* Block,
	   TStnEvent*     Event,
	   Int_t          (*F)(TStnDataBlock *, TStnEvent *, Int_t) = 0);

  virtual ~TStnNode();
					// ****** init methods

					// ****** accessors

  TStnDataBlock*   GetDataBlock       () { return  fObject; }
  TStnDataBlock**  GetDataBlockAddress() { return &fObject; }
  TBranch**        GetBranchAddress   () { return &fBranch; }
  TBranch*         GetBranch          () { return fBranch;  }
  TStnEvent*       GetEvent           () { return fEvent;   }

  virtual Int_t    GetEntry(Int_t Ientry);

					// ****** setters

  void             SetBranch(TBranch*   b    ) { fBranch = b;     }

					// ****** overloaded functions of 
					// TObject

  void Print(Option_t* option = "") const;

  ClassDef(TStnNode,0)
};

#endif
