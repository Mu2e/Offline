//-----------------------------------------------------------------------------
// Jan 10 2001 P.Murat: an I/O intelligent array of integers intended for
// fast I/O operations
// "I/O intelligence" means that when reading the data in the array, when it is 
// necessary, reallocates its data buffer
//_____________________________________________________________________________
#include "TBuffer.h"

#include "base/TStnArrayI.hh"


ClassImp (TStnArrayI)

//_____________________________________________________________________________
TStnArrayI::TStnArrayI() {

  fArray = NULL;
  fNBufferWords = 0;
  fNDataWords = 0;

}

//_____________________________________________________________________________
TStnArrayI::TStnArrayI(Int_t n) {

  fArray = NULL;
  fNBufferWords = 0;
  fNDataWords = 0;
  Set(n);

}

//_____________________________________________________________________________
TStnArrayI::TStnArrayI(Int_t* array, Int_t n) {

  fArray = NULL;
  fNBufferWords = 0;
  fNDataWords = 0;
  Set(array, n);

}

//_____________________________________________________________________________
TStnArrayI::TStnArrayI(TStnArrayI& rhs) {

  fArray = NULL;
  fNBufferWords = 0;
  fNDataWords = 0;
  Set(rhs.GetArray(), rhs.GetSize());

}

//_____________________________________________________________________________
TStnArrayI::~TStnArrayI() {
  delete[] fArray;
}

//_____________________________________________________________________________
Int_t TStnArrayI::Append(Int_t* array, Int_t nw) 
{
  Reserve(fNDataWords+nw);
  memcpy(fArray+fNDataWords,array,nw*sizeof(Int_t));
  fNDataWords += nw;
  return 0;
}

//_____________________________________________________________________________
Int_t TStnArrayI::Set(Int_t n) {
  
  Reserve(n);
  fNDataWords = n;

  return 0;
}

//_____________________________________________________________________________
Int_t TStnArrayI::Set(Int_t* array, Int_t nw) 
{
  Clear();
  return Append(array, nw);
}



//_____________________________________________________________________________
Int_t TStnArrayI::Reserve(Int_t nw) {
  
  if (fNBufferWords >= nw) return 0;
  return SetPhysicalSize( nw );

}

//_____________________________________________________________________________
Int_t TStnArrayI::SetPhysicalSize(Int_t nw) {

  // if logical size is larger, truncate it
  if(nw<fNDataWords) fNDataWords = nw;

  // expand in jumps
  if(nw>100) nw = (nw/100+1)*100;
  if(nw>1000) nw = (nw/1000+1)*1000;

  int* ia = new int[nw];

  // preserve the contents of the buffer  
  if (fArray) {
    memcpy(ia,fArray,fNDataWords*sizeof(Int_t));
    delete [] fArray;
  }
  fNBufferWords = nw;
  fArray = ia;

  return 0;
}

//_____________________________________________________________________________
void TStnArrayI::Print(const char* Opt) const {

  int loc = 0, ncolumns = 10;

  if (strcmp(Opt,"") != 0) {
    printf("----------------------------------------------------------------\n");
    printf(" %s\n",Opt);
    printf("----------------------------------------------------------------\n");
  }

  for (int i=0; i<fNDataWords; i++) {
    printf(" %12i",fArray[i]);
    loc++;
    if (loc == ncolumns) {
      printf("\n");
      loc = 0;
    }
  }
  if (loc != 0) printf("\n");
}

//_____________________________________________________________________________
void TStnArrayI::Streamer(TBuffer &R__b) {
  // Stream an object of class TStnArrayI. 

  Int_t n;
  if (R__b.IsReading()) {
    R__b >> n;
    if (n > 0) {
      Set( n );
      R__b.ReadFastArray(fArray,n);
    }
    else {
      Clear();
    }
  }   
  else {
    R__b << fNDataWords;
    if (fNDataWords > 0) {
      R__b.WriteFastArray(fArray,fNDataWords);
    }
  }
}


