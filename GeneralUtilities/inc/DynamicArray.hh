#ifndef GeneralUtilities_DynamicArray_hh
#define GeneralUtilities_DynamicArray_hh

//
//
// forms a dynamic array that you can size and resize
// 
//
// $Id: DynamicArray.hh,v 1.2 2011/05/17 15:41:35 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:41:35 $
//

#include <iostream>
#include <vector>
using namespace std;

template <typename T>
class DynamicArray
{
public:
  DynamicArray(){};

  DynamicArray(int rows, int cols): dArray(rows, vector<T>(cols)){}
  
  vector<T> & operator[](int i) 
  { 
    return dArray[i];
  }
 
  
  const vector<T> & operator[] (int i) const 
  { 
    return dArray[i];
  }
  
  void resize(int rows, int cols)//resize the two dimentional array .


  {
    dArray.resize(rows);
    for(int i = 0;i < rows;++i) dArray[i].resize(cols);
  }
private:
  vector<vector<T> > dArray;  
};
#endif /* GeneralUtilities_DynamicArray_hh */
