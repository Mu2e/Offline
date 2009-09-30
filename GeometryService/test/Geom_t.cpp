//
// Test features of the Geometry Service.
// 
// $Id: Geom_t.cpp,v 1.1 2009/09/30 22:57:47 kutschke Exp $
// $Author: kutschke $
// $Date: 2009/09/30 22:57:47 $
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
  catch (cms::Exception& x)
    {
      std::cerr << "cms::Exception caught:\n" << x << '\n';
      rc = 1;
    }
  catch (...)
    {
      std::cerr << "Unknown exception caught\n";
      rc = 2;
    }
  return rc;
}


