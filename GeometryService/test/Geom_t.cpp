//
// Test features of the Geometry Service.
// 
// $Id: Geom_t.cpp,v 1.2 2011/05/17 15:36:00 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:36:00 $
//
// Original author Jim Kowalkowski.
//
#include <cassert>
#include <iostream>
#include <string>
#include <vector>

#include "GeometryService/inc/GeometryService.hh"

int main()
{
  int rc = 0;
  try
    {
      // rc = work();
    }
  catch (cet::exception& x)
    {
      std::cerr << "cet::exception caught:\n" << x << '\n';
      rc = 1;
    }
  catch (...)
    {
      std::cerr << "Unknown exception caught\n";
      rc = 2;
    }
  return rc;
}


