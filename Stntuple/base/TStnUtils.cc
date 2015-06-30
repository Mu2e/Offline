///////////////////////////////////////////////////////////////////////////////
// 2003-12-30 P.Murat
///////////////////////////////////////////////////////////////////////////////

#include <cstdlib>
#include "Stntuple/base/TStnUtils.hh"


ClassImp(TStnUtils)

//_____________________________________________________________________________
TStnUtils::TStnUtils() {
}

//_____________________________________________________________________________
TStnUtils::~TStnUtils() {
}


//_____________________________________________________________________________
Int_t TStnUtils::ReadTable(const char* InputFile , 
			   TObjArray*  Table     , 
			   const char* Delimitors) 
{
  // assume table has enough length to hold all the columns

  TStnArrayF* column;

  char c[200], *token;

  float x[1];

  int ind, line = 0, ncolumns = 0, nmax=200;

  FILE* file = fopen(InputFile,"r");
  while ((c[0]=getc(file)) != EOF) {

					// check if it is a comment line
    ind = 0;
    if (c[0] != '#') {
      ungetc(c[0],file);
      fgets(c,nmax,file);
//-----------------------------------------------------------------------------
// parse single line
//-----------------------------------------------------------------------------
      //      printf("%s",c);
      token = strtok (c,Delimitors);
      do {
//  	printf("new token:%s\n",token);
//  	if (getchar() == 'q') return;
	x[0] = atof(token);
	if (line == 0) {	// allocate vector for each column

	  column = new TStnArrayF(100);
	  Table->Add((TObject*) column);
	  ncolumns++;
	}
	else {			// number of columns is already known

	  column = (TStnArrayF*) Table->UncheckedAt(ind);
	}
	column->Append(x,1);
	ind++;
      } while ((token = strtok (0,Delimitors)));

      line++;
    }
    else {
				// skip line
      fgets(c,nmax,file);
    }
  }
  return 0;
}

//-----------------------------------------------------------------------------
// read array: assume ASCII file formatted arbitrary way, but read all the
// numbers (separated by " ,") - integers - into one array
// skip lines starting from "#"
//-----------------------------------------------------------------------------
Int_t TStnUtils::ReadArrayI(const char* InputFile , 
			    TStnArrayI* Array     , 
			    const char* Delimitors) 
{
  char c[200], *token;

  int   x[1];

  int nmax=200;

  FILE* file = fopen(InputFile,"r");
  while ((c[0]=getc(file)) != EOF) {
    if (c[0] != '#') {
      ungetc(c[0],file);
      fgets(c,nmax,file);
//-----------------------------------------------------------------------------
// parse single line
//-----------------------------------------------------------------------------
      //      printf("%s",c);
      token = strtok (c,Delimitors);
      do {

//  	printf("new token:%s\n",token);
//  	if (getchar() == 'q') return -1;

	x[0] = atoi(token);
	Array->Append(x,1);
      } while ((token = strtok (0,Delimitors)));
    }
    else {
				// skip line
      fgets(c,nmax,file);
    }
  }

  return 0;
}

//-----------------------------------------------------------------------------
// read array: assume ASCII file formatted arbitrary way, but read all the
// numbers (separated by " ,") - integers - into one array
// skip lines starting from "#"
//-----------------------------------------------------------------------------
Int_t TStnUtils::ReadArrayF(const char* InputFile , 
			    TStnArrayF* Array     , 
			    const char* Delimitors) 
{
  char c[200], *token;

  float  x[1];

  int nmax=200;

  FILE* file = fopen(InputFile,"r");
  while ((c[0]=getc(file)) != EOF) {
    if (c[0] != '#') {
      ungetc(c[0],file);
      fgets(c,nmax,file);
//-----------------------------------------------------------------------------
// parse single line
//-----------------------------------------------------------------------------
      //      printf("%s",c);
      token = strtok (c,Delimitors);
      do {

// 		printf("new token:%s\n",token);
// 		if (getchar() == 'q') return;

	x[0] = atof(token);
	Array->Append(x,1);
      } while ((token = strtok (0,Delimitors)));
    }
    else {
				// skip line
      fgets(c,nmax,file);
    }
  }

  return 0;
}

