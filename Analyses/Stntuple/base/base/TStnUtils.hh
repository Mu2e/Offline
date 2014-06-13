///////////////////////////////////////////////////////////////////////////////
// Dec 30 2003 P.Murat
///////////////////////////////////////////////////////////////////////////////
#include "TObject.h"
#include "TObjArray.h"

#include "TStnArrayI.hh"
#include "TStnArrayF.hh"

class TStnUtils: public TObject {
public:
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//  functions:
//-----------------------------------------------------------------------------
public:
  TStnUtils();
  ~TStnUtils();
//-----------------------------------------------------------------------------
// other methods: 
// table: N columns of floats
//-----------------------------------------------------------------------------
  static   int ReadTable (const char* InputFile, 
			  TObjArray*  Table    , 
			  const char* Delimitors = " ,");

  static   int ReadArrayI(const char* InputFile, 
			  TStnArrayI* Array    ,
			  const char* Delimitors = " ,");

  static   int ReadArrayF(const char* InputFile, 
			  TStnArrayF* Array    ,
			  const char* Delimitors = " ,");
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  ClassDef(TStnUtils,0)
};
