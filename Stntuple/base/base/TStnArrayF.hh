//-----------------------------------------------------------------------------
//  Jan 10 2001 P.Murat: I/O intelligent arary of integers
//-----------------------------------------------------------------------------
#ifndef TStnArrayF_hh
#define TStnArrayF_hh

#include <iostream>
#include "TNamed.h"

class TStnArrayF: public TNamed {
  Int_t      fNDataWords;		// # of real data words in the buffer
  Int_t      fNBufferWords;		// # size of the buffer
  Float_t*   fArray;		        //!# array of data

public:
					// ****** constructors and destructor
  TStnArrayF();
  TStnArrayF(Int_t n);   // set physical and logical size to n
  TStnArrayF(Float_t* array, Int_t n);   // set phys. and log. size to n, fill
  TStnArrayF(TStnArrayF& rhs);
  virtual ~TStnArrayF();

  Int_t         NDataWords()       { return fNDataWords; }
  Int_t         NDataWords() const { return fNDataWords; }
  Int_t         GetSize()    const { return fNDataWords; }
  Int_t         BufferSize()       { return fNBufferWords; }
  Float_t*        GetArray()         { return fArray; }
  const Float_t*  GetArray() const   { return fArray; }

  // append to logical contents (physical size will be increased as necessary)
  Int_t         Append(Float_t* array, int nwords);

  // set logical size (physical size will be increased as necessary)
  Int_t         Set(Int_t n);
  // set logical size and set contents
  Int_t         Set(Float_t* array, Int_t n);

  // using these does not change physical or logical size
  // references outside logical size gives errors
  Float_t       At(Int_t i) const;
  Float_t       operator[](Int_t i) const;
  Float_t      &operator[](Int_t i);

  // set physical size of the array to max(n,current size), 
  // saves current contents
  Int_t         Reserve(Int_t n);

  // set physical size of the array to n, 
  // may truncate logical size
  Int_t         SetPhysicalSize(Int_t n);

  void          Clear(Option_t* opt = "") { fNDataWords = 0; }
  void          Print(Option_t* opt = "") const ;

  ClassDef(TStnArrayF,2)		// an intelligent array of integers
};



inline Float_t TStnArrayF::At(Int_t i) const
{
  if(i<0 || i>=fNDataWords){
    TObject::Error("operator[]",
	"access of entry %d, while limits are [%d,%d]\n",i,0,fNDataWords-1);
  }
   return fArray[i];
}

inline Float_t &TStnArrayF::operator[](Int_t i) 
{
  if(i<0 || i>=fNDataWords){
    TObject::Error("operator[]",
	"access of entry %d, while limits are [%d,%d]\n",i,0,fNDataWords-1);
    i = 0;
  }
  return fArray[i];
}

inline Float_t TStnArrayF::operator[](Int_t i) const
{
  if(i<0 || i>=fNDataWords){
    TObject::Error("operator[]",
	"access of entry %d, while limits are [%d,%d]\n",i,0,fNDataWords-1);
    return 0;
  }
  return fArray[i];
}


#endif

