//
//     c++ rewrite of the Karimaki circle fit (CLEFIT)
//     original code CERN acbz.f
//
// $Id: KarimakiCircle.cc,v 1.1 2012/05/15 07:51:36 tassiell Exp $
// $Author: tassiell $
// $Date: 2012/05/15 07:51:36 $
//
// Original author G. Tassielli
//

// Mu2e includes.
#include "FastPatternReco/inc/KarimakiCircle.hh"

using namespace std;

namespace mu2e {

/* **************************************************************** */
int trasa3_(float *a, float *b)
{
        float c_[6];
        long int i_;
        float e1, e2, e3, e4, e5, e6, e7, e8, e9;

        /* ------------------------------------------------------------- * */
        /*     MATRIX OPERATION ABA'-->B FOR 3X3 MATRICES                * */
        /*     B is a packed symmetric matrix                            * */
        /*     A is 3x3 matrix packed row-wise     (A' means transpose)  * */
        /* ------------------------------------------------------------- * */
        /* Parameter adjustments */
        --b;
        --a;

        /* Function Body */
        for (i_ = 1; i_ <= 6; ++i_) {
                c_[i_ - 1] = b[i_];
        }
        e1 = a[1] * c_[0] + a[2] * c_[1] + a[3] * c_[3];
        e2 = a[1] * c_[1] + a[2] * c_[2] + a[3] * c_[4];
        e3 = a[1] * c_[3] + a[2] * c_[4] + a[3] * c_[5];
        e4 = a[4] * c_[0] + a[5] * c_[1] + a[6] * c_[3];
        e5 = a[4] * c_[1] + a[5] * c_[2] + a[6] * c_[4];
        e6 = a[4] * c_[3] + a[5] * c_[4] + a[6] * c_[5];
        e7 = a[7] * c_[0] + a[8] * c_[1] + a[9] * c_[3];
        e8 = a[7] * c_[1] + a[8] * c_[2] + a[9] * c_[4];
        e9 = a[7] * c_[3] + a[8] * c_[4] + a[9] * c_[5];
        b[1] = a[1] * e1 + a[2] * e2 + a[3] * e3;
        b[2] = a[4] * e1 + a[5] * e2 + a[6] * e3;
        b[3] = a[4] * e4 + a[5] * e5 + a[6] * e6;
        b[4] = a[7] * e1 + a[8] * e2 + a[9] * e3;
        b[5] = a[7] * e4 + a[8] * e5 + a[9] * e6;
        b[6] = a[7] * e7 + a[8] * e8 + a[9] * e9;
        return 0;
} /* trasa3_ */
/* **************************************************************** */
KarimakiCircle::KarimakiCircle() : rho(0.0), phi(0.0), dca(0.0), chicir(0.0), xpca(0.0), ypca(0.0), xx0(0.0), yy0(0.0), cosf(0.0), sinf(0.0), dirtx(0.0), dirty(0.0) {

        covrfd[0] = covrfd[1] = covrfd[2] = covrfd[3] = covrfd[4] = covrfd[5] = 0.0;

        ierror = 0;
        s1 = 0.0;
        s2 = 0.0;
        s3 = 0.0;
        s4 = 0.0;
        s5 = 0.0;
        s6 = 0.0;
        s7 = 0.0;
        s8 = 0.0;
        s9 = 0.0;
}

void KarimakiCircle::addHit ( float xx, float yy, float errxx, float erryy ){
        points.push_back( circPoint(xx, yy, errxx, erryy) );
        if ( points.size()==1 ){
                xx0 = points.back()._xx;
                yy0 = points.back()._yy;
        }

        /* --- calculate sums for fit */
        xc = points.back()._xx - xx0;
        yc = points.back()._yy - yy0;
        wt = 1.0/( points.back()._errxx * points.back()._erryy );
        wx = wt * xc;
        wy = wt * yc;
        rr = xc * xc + yc * yc;
        wr = wt * rr;
        s1 += wt;
        s2 += wx;
        s3 += wy;
        s4 += wx * xc;
        s5 += wx * yc;
        s6 += wy * yc;
        s7 += wx * rr;
        s8 += wy * rr;
        s9 += wr * rr;
        if (points.size()<3 || s1 <= 0.0) {
                ierror = 1;
        }
        ierror = 0;
}

void KarimakiCircle::testHit ( float xx, float yy, float errxx, float erryy ){

        /* --- calculate sums for fit */
        xc = xx - xx0;
        yc = yy - yy0;
        wt = 1.0/( errxx * erryy );
        wx = wt * xc;
        wy = wt * yc;
        rr = xc * xc + yc * yc;
        wr = wt * rr;
        s1 += wt;
        s2 += wx;
        s3 += wy;
        s4 += wx * xc;
        s5 += wx * yc;
        s6 += wy * yc;
        s7 += wx * rr;
        s8 += wy * rr;
        s9 += wr * rr;
        ierror = 0;
}

bool KarimakiCircle::computeBestCirc( bool useOldDir, float pointDirX, float pointDirY, float txref, float tyref ) {
        if (ierror!=0) return false;

        /* --- Solve the fitted parameters */

        s1i = 1.0 / s1;
        sr = s4 + s6;
        hsr = sr * 0.5;
        xmean = s1i * s2;
        ymean = s1i * s3;
        rrmean = s1i * sr;
        cov1 = s1i * (s4 - s2 * xmean);
        cov2 = s1i * (s5 - s2 * ymean);
        cov3 = s1i * (s6 - s3 * ymean);
        cov4 = s1i * (s7 - s2 * rrmean);
        cov5 = s1i * (s8 - s3 * rrmean);
        cov6 = s1i * (s9 - sr * rrmean);
        if (cov6 <= 0.0) {
                ierror = 1;
                return false;

        }
        y2fi = (cov2 * cov6 - cov4 * cov5) * 2.0;
        x2fi = cov6 * (cov1 - cov3) - cov4 * cov4 + cov5 * cov5;
        fifit = std::atan2(y2fi, x2fi) * 0.5;
        cosf = std::cos(fifit);
        sinf = std::sin(fifit);
        hapfit = (sinf * cov4 - cosf * cov5) / cov6;
        delfit = -hapfit * rrmean + sinf * xmean - cosf * ymean;
        apu = std::sqrt( 1.0 - hapfit * 4.0 * delfit);
        rhof = hapfit * 2.0 / apu;
        dft = delfit * 2.0 / (apu + 1.);
        rod1 = 1.0 / apu;
        rod2 = rod1 * rod1;
        sinf2 = sinf * sinf;
        cosf2 = cosf * cosf;
        sinff = sinf * 2.0 * cosf;
        sa = sinf * s2 - cosf * s3;
        saa = sinf2 * s4 - sinff * s5 + cosf2 * s6;
        sxyr = sinf * s7 - cosf * s8;

        rho = rhof;
        phi = fifit;
        dca = dft;
        chicir = rod2 * (-delfit * sa - hapfit * sxyr + saa);
        xpca = xx0 + dca * sinf;
        ypca = yy0 - dca * cosf;

        ierror = 0;

        /* --- Error estimation ro,fi,d */

        //if (*mode > 0) {
        sb = cosf * s2 + sinf * s3;
        sg = (sinf2 - cosf2) * s5 + sinf * cosf * (s4 - s6);
        w1 = s9 * 0.25 - dft * (sxyr - dft * (saa + hsr - dft *
                        (sa - dft * 0.25 * s1)));
        w2 = -rod1 * ((cosf * s7 + sinf * s8) * 0.5 -
                        dft * (sg - dft * 0.5 * sb));
        w3 = rod2 * (cosf2 * s4 + sinff * s5 + sinf2 * s6);
        w4 = rhof * (sxyr * (-0.5) + dft * saa) + rod1 * hsr - dft * 0.5 *
                        ((rod1 * 2.0 + rhof * dft) * sa - dft * rod1 * s1);
        w5 = rod1 * rhof * sg - rod2 * sb;
        w6 = rhof * (rhof * saa - rod1 * 2.0 * sa) + rod2 * s1;
        sd1 = w3 * w6 - w5 * w5;
        sd2 = -w2 * w6 + w4 * w5;
        sd3 = w2 * w5 - w3 * w4;
        detinv = 1.0 / (w1 * sd1 + w2 * sd2 + w4 * sd3);
        covrfd[0] = detinv * sd1;
        covrfd[1] = detinv * sd2;
        covrfd[2] = detinv * (w1 * w6 - w4 * w4);
        covrfd[3] = detinv * sd3;
        covrfd[4] = detinv * (w2 * w4 - w1 * w5);
        covrfd[5] = detinv * (w1 * w3 - w2 * w2);
        xdero = (rhof * sxyr - rod1 * 2.0 * saa + (rod1 + 1.0) * dft * sa) * 0.5;
        edero = dft * xdero;
        ededi = rhof * xdero;
        drof = covrfd[0] * edero + covrfd[3] * ededi;
        dfif = covrfd[1] * edero + covrfd[4] * ededi;
        ddif = covrfd[3] * edero + covrfd[5] * ededi;
        rho += drof;
        dca += ddif;
        phi += dfif;
        /* Computing 2nd power */
        r_1 = rho * dca + 1.0;
        chicir = r_1 * r_1 * chicir / rod2;
        //}

        if(useOldDir) {
                PropagateToRef( txref, tyref, false );
        } else {
                if (pointDirX!=0.0 && pointDirY!=0.0) {
                        dirtx = xx0 - pointDirX; //circfi1.xx0 - xx[1];
                        dirty = yy0 - pointDirY; //circfi1.yy0 - yy[1];
                }
                else {
                        dirtx = xx0 - points.at(1)._xx;//circfi1.xx0 - xx[1];
                        dirty = yy0 - points.at(1)._yy;//circfi1.yy0 - yy[1];
                        //int m3 = np()/3+1;
                        //dirtx = points.at(m3)._xx - xx0; //circfi1.xx0 - xx[1];
                        //dirty = points.at(m3)._yy - yy0; //circfi1.yy0 - yy[1];
                }

                PropagateToRef( txref, tyref, false );
        }
        return true;
}

void KarimakiCircle::PropagateToRef(float &txref, float &tyref, bool notFirst ) {

        if (notFirst) {
                dirtx = std::cos(phi);
                dirty = std::sin(phi);
                sinf = std::sin(phi);
                cosf = std::cos(phi);
        }

        /* --- Propagate parameters to  XREF,YREF */

        xref  = txref;
        yref  = tyref;
        xmove = xpca - xref;
        ymove = ypca - yref;
        rod1 = rho * dca + 1.0;
        dperp = xmove * sinf - ymove * cosf;
        dpara = xmove * cosf + ymove * sinf;
        zee = dperp * dperp + dpara * dpara;
        aa = dperp * 2.0 + rho * zee;
        uu = std::sqrt(rho * aa + 1.0);
        sq1ai = 1.0 / (uu + 1.0);
        bb = rho * xmove + sinf;
        cc = -rho * ymove + cosf;
        phi = std::atan2(bb, cc);
        dca = aa * sq1ai;

        /* --- Propagate error matrix to XREF,YREF */

        //if (*mode > 0) {
        vv = rho * dperp + 1.0;
        xee = 1.0 / (cc * cc + bb * bb);
        xla = aa * 0.5 * (sq1ai * sq1ai) / uu;
        xmu = sq1ai / uu + rho * xla;
        //xjacob[9];
        xjacob[0] = 1.0;
        xjacob[1] = 0.0;
        xjacob[2] = 0.0;
        xjacob[3] = xee * dpara;
        xjacob[4] = xee * rod1 * vv;
        xjacob[5] = -xjacob[3] * (rho * rho);
        xjacob[6] = xmu * zee - xla * aa;
        xjacob[7] = xmu * 2.0 * rod1 * dpara;
        xjacob[8] = xmu * 2.0 * vv;
        trasa3_(xjacob, covrfd);
        sinf = std::sin(phi);
        cosf = std::cos(phi);
        //}

        /* --- check direction */

        dirtes = cosf * dirtx + sinf * dirty;
        if (dirtes < 0.0) {
                phi += picons;
                cosf = -cosf;
                sinf = -sinf;
                dca = -dca;
                rho = -rho;
                //if (imod > 0) {
                covrfd[1] = -covrfd[1];
                covrfd[4] = -covrfd[4];
                //}
        }
        if (phi >= pi2con) {
                phi -= pi2con;
        }
        if (phi < 0.0) {
                phi += pi2con;
        }
        xpca = xref + dca * sinf;
        ypca = yref - dca * cosf;

}

bool KarimakiCircle::summCirc( KarimakiCircle const & addKrm, int skipFrstNPnts ) {
        cout<<"Krm::Npoints before sum "<<points.size()<<endl;
        for ( std::vector<circPoint>::const_iterator addPoints_it = (addKrm.points.begin()+skipFrstNPnts); addPoints_it != addKrm.points.end(); ++addPoints_it) {
                addHit(addPoints_it->_xx, addPoints_it->_yy, addPoints_it->_errxx, addPoints_it->_erryy);
        }
        cout<<"Krm::Npoints after sum "<<points.size()<<endl;
        return computeBestCirc(false,0.0,0.0,xref,yref);
}

bool KarimakiCircle::bestmergeCirc( KarimakiCircle const & addKrm, int skipFrstNPnts ) {
        cout<<"Krm::Npoints before sum "<<points.size()<<endl;
        float ndof = points.size() -3;
        //float startChi2 = chicir/ndof;
        float tmpChi2 = 0.0;
        for ( std::vector<circPoint>::const_iterator addPoints_it = (addKrm.points.begin()+skipFrstNPnts); addPoints_it != addKrm.points.end(); ++addPoints_it) {
                testHit( addPoints_it->_xx, addPoints_it->_yy, addPoints_it->_errxx, addPoints_it->_erryy );
                if ( computeBestCirc(false,0.0,0.0,xref,yref) ) {
                        ++ndof;
                        tmpChi2 = chicir/ndof;
                        if ( tmpChi2>0.001 && tmpChi2<10.0 ) {
                                points.push_back( *addPoints_it );
                        } else {
                                --ndof;
                                testHit( addPoints_it->_xx, addPoints_it->_yy, -addPoints_it->_errxx, addPoints_it->_erryy );
                        }
                }
                else {
                        testHit( addPoints_it->_xx, addPoints_it->_yy, -addPoints_it->_errxx, addPoints_it->_erryy );
                }
        }
        cout<<"Krm::Npoints after sum "<<points.size()<<endl;
        return computeBestCirc(false,0.0,0.0,xref,yref);
}

bool KarimakiCircle::mergeHitOfCirc( KarimakiCircle const & addKrm, size_t hitToAdd ) {
        float ndof = points.size() -3;
        float startDistChi2 = fabs(1.0-chicir/ndof);
        float tmpChi2 = 0.0, tmpDistChi2;
        circPoint const& addPoints = addKrm.points.at(hitToAdd);
        bool added = false;
        testHit( addPoints._xx, addPoints._yy, addPoints._errxx, addPoints._erryy );
        if ( computeBestCirc(false,0.0,0.0,xref,yref) ) {
                ++ndof;
                tmpChi2 = chicir/ndof;
                tmpDistChi2 = fabs(1.0-tmpChi2);
                if ( tmpChi2>0.001 && tmpChi2<10.0 && tmpDistChi2<startDistChi2 ) {
                        points.push_back( addPoints );
                        startDistChi2=tmpDistChi2;
                        added = true;
                } else {
                        --ndof;
                        testHit( addPoints._xx, addPoints._yy, -addPoints._errxx, addPoints._erryy );
                }
        }
        else {
                testHit( addPoints._xx, addPoints._yy, -addPoints._errxx, addPoints._erryy );
        }
        computeBestCirc(false,0.0,0.0,xref,yref);
        return added;
}

void KarimakiCircle::removeHit (size_t iHit) {
        if (iHit==0) {
                points.erase(points.begin());
                covrfd[0] = covrfd[1] = covrfd[2] = covrfd[3] = covrfd[4] = covrfd[5] = 0.0;

                ierror = 0;
                s1 = 0.0;
                s2 = 0.0;
                s3 = 0.0;
                s4 = 0.0;
                s5 = 0.0;
                s6 = 0.0;
                s7 = 0.0;
                s8 = 0.0;
                s9 = 0.0;
                xx0 = points.at(0)._xx;
                yy0 = points.at(0)._yy;

                /* --- calculate sums for fit */
                for ( std::vector<circPoint>::iterator points_it = (points.begin()+1); points_it != points.end(); ++points_it ) {
                        xc = points_it->_xx - xx0;
                        yc = points_it->_yy - yy0;
                        wt = 1.0/( points_it->_errxx * points_it->_erryy );
                        wx = wt * xc;
                        wy = wt * yc;
                        rr = xc * xc + yc * yc;
                        wr = wt * rr;
                        s1 += wt;
                        s2 += wx;
                        s3 += wy;
                        s4 += wx * xc;
                        s5 += wx * yc;
                        s6 += wy * yc;
                        s7 += wx * rr;
                        s8 += wy * rr;
                        s9 += wr * rr;
                }
                computeBestCirc(true);
        } else {
                std::vector<circPoint>::iterator removeKrmPnt_it = points.begin()+iHit;
                testHit( removeKrmPnt_it->_xx, removeKrmPnt_it->_yy, -removeKrmPnt_it->_errxx, removeKrmPnt_it->_erryy);
                points.erase(removeKrmPnt_it);
                computeBestCirc(true,0.0,0.0,xref,yref);
        }

}

KarimakiCircle & KarimakiCircle::operator += ( KarimakiCircle const & addKrm ) {
        cout<<"Krm::Npoints before sum "<<points.size()<<endl;
        for ( std::vector<circPoint>::const_iterator addPoints_it = addKrm.points.begin(); addPoints_it != addKrm.points.end(); ++addPoints_it) {
                addHit(addPoints_it->_xx, addPoints_it->_yy, addPoints_it->_errxx, addPoints_it->_erryy);
        }
        cout<<"Krm::Npoints after sum "<<points.size()<<endl;
        computeBestCirc(false,0.0,0.0,xref,yref);

        return *this;
}

KarimakiCircle & KarimakiCircle::operator= ( KarimakiCircle const& tmpKrm ) {
        points = tmpKrm.points;
        rho = tmpKrm.rho;
        phi = tmpKrm.phi;
        dca = tmpKrm.dca;
        chicir = tmpKrm.chicir;
        xpca = tmpKrm.xpca;
        ypca = tmpKrm.ypca;
        covrfd[0] = tmpKrm.covrfd[0];
        covrfd[1] = tmpKrm.covrfd[1];
        covrfd[2] = tmpKrm.covrfd[2];
        covrfd[3] = tmpKrm.covrfd[3];
        covrfd[4] = tmpKrm.covrfd[4];
        covrfd[5] = tmpKrm.covrfd[5];
        xx0 = tmpKrm.xx0;
        yy0 = tmpKrm.yy0;

        xc = tmpKrm.xc;
        yc = tmpKrm.yc;
        wt = tmpKrm.wt;
        wx = tmpKrm.wx;
        wy = tmpKrm.wy;
        rr = tmpKrm.rr;
        wr = tmpKrm.wr;
        s1 = tmpKrm.s1;
        s2 = tmpKrm.s2;
        s3 = tmpKrm.s3;
        s4 = tmpKrm.s4;
        s5 = tmpKrm.s5;
        s6 = tmpKrm.s6;
        s7 = tmpKrm.s7;
        s8 = tmpKrm.s8;
        s9 = tmpKrm.s9;
        ierror = tmpKrm.ierror;

        xref = tmpKrm.xref;
        yref = tmpKrm.yref;

}

}  // end namespace mu2e
