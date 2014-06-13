//-----------------------------------------------------------------------------
//  Jan 10 2001 P.Murat: I/O intelligent arary of integers
//-----------------------------------------------------------------------------
#ifndef TStnArrayI_hh
#define TStnArrayI_hh

#include <iostream>
#include "TNamed.h"

class TStnArrayI: public TNamed {
  Int_t      fNDataWords;		// # of real data words in the buffer
  Int_t      fNBufferWords;		// # size of the buffer
  Int_t*     fArray;		        //!# array of data

public:
					// ****** constructors and destructor
  TStnArrayI();
  TStnArrayI(Int_t n);   // set physical and logical size to n
  TStnArrayI(Int_t* array, Int_t n);   // set phys. and log. size to n, fill
  TStnArrayI(TStnArrayI& rhs);
  virtual ~TStnArrayI();

  Int_t         NDataWords()       { return fNDataWords; }
  Int_t         NDataWords() const { return fNDataWords; }
  Int_t         GetSize()    const { return fNDataWords; }
  Int_t         BufferSize()       { return fNBufferWords; }
  Int_t*        GetArray()         { return fArray; }
  const Int_t*  GetArray() const   { return fArray; }

  // append to logical contents (physical size will be increased as necessary)
  Int_t         Append(Int_t* array, int nwords);

  // set logical size (physical size will be increased as necessary)
  Int_t         Set(Int_t n);
  // set logical size and set contents
  Int_t         Set(Int_t* array, Int_t n);

  // using these does not change physical or logical size
  // references outside logical size gives errors
  Int_t         At(Int_t i) const;
  Int_t         operator[](Int_t i) const;
  Int_t        &operator[](Int_t i);

  // set physical size of the array to max(n,current size), 
  // saves current contents
  Int_t         Reserve(Int_t n);

  // set physical size of the array to n, 
  // may truncate logical size
  Int_t         SetPhysicalSize(Int_t n);

  void          Clear(Option_t* opt = "") { fNDataWords = 0; }
  void          Print(Option_t* opt = "") const ;

  ClassDef(TStnArrayI,2)		// an intelligent array of integers
};



inline Int_t TStnArrayI::At(Int_t i) const
{
  if(i<0 || i>=fNDataWords){
    TObject::Error("operator[]",
	"access of entry %d, while limits are [%d,%d]\n",i,0,fNDataWords-1);
  }
   return fArray[i];
}

inline Int_t &TStnArrayI::operator[](Int_t i) 
{
  if(i<0 || i>=fNDataWords){
    TObject::Error("operator[]",
	"access of entry %d, while limits are [%d,%d]\n",i,0,fNDataWords-1);
    i = 0;
  }
  return fArray[i];
}

inline Int_t TStnArrayI::operator[](Int_t i) const
{
  if(i<0 || i>=fNDataWords){
    TObject::Error("operator[]",
	"access of entry %d, while limits are [%d,%d]\n",i,0,fNDataWords-1);
    return 0;
  }
  return fArray[i];
}


#endif

