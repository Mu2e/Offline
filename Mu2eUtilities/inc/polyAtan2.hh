#ifndef Mu2eUtilities_polyAtan2_hh
#define Mu2eUtilities_polyAtan2_hh
//
//
//
// Original author Giani Pezzullo
//

// Modern CLHEP
#include "CLHEP/Vector/ThreeVector.h"


namespace mu2e {

  inline const float     ApproxAtan(float z)
  {
    const float n1 = 0.97239411f;
    const float n2 = -0.19194795f;
    return (n1 + n2 * z * z) * z;
  }

  inline const float     polyAtan2(float y, float x){
    float   PI_2 = M_PI/2.;
    if (x != 0.0f)
      {
        if (fabsf(x) > fabsf(y))
          {
            const float z = y / x;
            if (x > 0.0)
              {
                // atan2(y,x) = atan(y/x) if x > 0
                return ApproxAtan(z);
              }
            else if (y >= 0.0)
              {
                // atan2(y,x) = atan(y/x) + PI if x < 0, y >= 0
                return ApproxAtan(z) + M_PI;
              }
            else
              {
                // atan2(y,x) = atan(y/x) - PI if x < 0, y < 0
                return ApproxAtan(z) - M_PI;
              }
          }
        else // Use property atan(y/x) = PI/2 - atan(x/y) if |y/x| > 1.
          {
            const float z = x / y;
            if (y > 0.0)
              {
                // atan2(y,x) = PI/2 - atan(x/y) if |y/x| > 1, y > 0
                return -ApproxAtan(z) + PI_2;
              }
            else
              {
                // atan2(y,x) = -PI/2 - atan(x/y) if |y/x| > 1, y < 0
                return -ApproxAtan(z) - PI_2;
              }
          }
      }
    else
      {
        if (y > 0.0f) // x = 0, y > 0
          {
            return PI_2;
          }
        else if (y < 0.0f) // x = 0, y < 0
          {
            return -PI_2;
          }
      }
    return 0.0f; // x,y = 0. Could return NaN instead.
  }

}
#endif /* Mu2eUtilities_polyAtan2_hh */
