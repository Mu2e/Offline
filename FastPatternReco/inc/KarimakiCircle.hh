#ifndef FastPatternReco_KarimakiCircle_hh
#define FastPatternReco_KarimakiCircle_hh
//
//     c++ rewrite of the Karimaki circle fit (CLEFIT)
//     original code CERN acbz.f
//
// $Id: KarimakiCircle.hh,v 1.2 2012/05/16 19:30:07 tassiell Exp $
// $Author: tassiell $
// $Date: 2012/05/16 19:30:07 $
//
// Original author G. Tassielli
//

/* ------------------------------------------------------------------ * */
/*     Non-iterative circle fit       (V. Karimaki/1991)              * */
/*                                                                    * */
/*     XREF,YREF = reference point coordinates           (IN)         * */
/*     XX,YY     = arrays of measured coordinates        (IN)         * */
/*     WW        = array of weigths                      (IN)         * */
/*     NP        = number of points                      (IN)         * */
/*     IERROR    = error flag (=0 if fit OK)            (OUT)         * */
/*                                                                    * */
/*     Fit results in COMMON/CIRCFI/:                                 * */
/*     RHO       = fitted curvature                     (OUT)         * */
/*     PHI       = fitted direction                     (OUT)         * */
/*     DCA       = fitted distance to (XREF,YREF)       (OUT)         * */
/*     CHICIR    = chi squared of the fit               (OUT)         * */
/*     XPCA,YPCA = point on circle closest to XREF,YREF (OUT)         * */
/*     COVRFD(6) = covariance matrix of RHO,PHI,DCA     (OUT)         * */
/*                 in lower triangular form                           * */
/*                                                                    * */
/*    NOTES:  1. Weights should be negative for point removal         * */
/* ------------------------------------------------------------------ * */

// C++ includes.
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <vector>

// CLHEP includes.
#include "CLHEP/Units/PhysicalConstants.h"

using namespace std;

namespace mu2e {

/* Table of constant values */
double picons = CLHEP::pi;//(float)3.14159265;
double pi2con = CLHEP::twopi;//(float)6.28318531;

struct circPoint {
        circPoint (float &xx, float &yy, float &errxx, float &erryy ) :
                _xx(xx),
                _yy(yy),
                _errxx(errxx),
                _erryy(erryy)
        {}
        float _xx;
        float _yy;
        float _errxx;
        float _erryy;
};

class KarimakiCircle {
public:
        KarimakiCircle();
        ~KarimakiCircle() {}

        void addHit ( float xx, float yy, float errxx, float erryy );
        bool checkBfrAddHit( float xx, float yy, float errxx, float erryy, float maxChi2=10.0, int minNDOFcutImprvChi2=20 );
        bool rejectHits(std::vector<circPoint>::iterator &points_it, float maxChi2=10.0);
        void testHit ( float xx, float yy, float errxx, float erryy );

        size_t np() const{
                return points.size();
        }

        bool computeBestCirc( bool useOldDir=false, float pointDirX=0.0, float pointDirY=0.0, float xref=0.0, float yref=0.0 );

        void PropagateToRef(float &xref, float &yref, bool notFirst=true );

        bool summCirc( KarimakiCircle const & addKrm, int skipFrstNPnts=0 );
        bool bestmergeCirc( KarimakiCircle const & addKrm, int skipFrstNPnts=0 );
        bool mergeHitOfCirc( KarimakiCircle const & addKrm, size_t hitToAdd );
        void removeHit (size_t iHit);
        KarimakiCircle & operator += ( KarimakiCircle const & addKrm );
        KarimakiCircle & operator= ( KarimakiCircle const& tmpKrm );

        std::vector<circPoint> points;
        float rho, phi, dca, chicir, xpca, ypca, covrfd[6], xx0, yy0;
        long int ierror;

private:

        double s1, s2, s3, s4, s5, s6, s7, s8, s9;
        float cosf, sinf;
        float dirtx, dirty, dirtes;

        float xc, yc, wt, wx, wy, rr, wr;
        float s1i, sr, hsr, xmean, ymean, rrmean;
        float cov1, cov2, cov3, cov4, cov5, cov6;
        float y2fi, x2fi, fifit, hapfit, delfit, apu;
        float rhof, dft, rod1, rod2, sinf2, cosf2, sinff;
        float sa, saa, sxyr, sb, sg, sd1, sd2, sd3;
        float w1, w2, w3, w4, w5, w6;
        float detinv;
        float xdero, edero, ededi, drof, dfif, ddif;
        float r_1;
        float xmove, ymove, dperp, dpara;
        float zee, aa, uu, sq1ai, bb, cc;
        float vv, xee, xla, xmu;
        float xjacob[9];
        float xref, yref;

};

}  // end namespace mu2e

#endif /* FastPatternReco_KarimakiCircle_hh */
