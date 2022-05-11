
// C++ includes.
#include <iostream>
#include <string>
#include <sstream>

#include "Offline/ParticleID/inc/PIDUtilities.hh"

using namespace std;

namespace mu2e {


  TH1D *PIDUtilities::th1dmorph(
                TH1D *hist1=0,TH1D *hist2=0,
                double par1=0,double par2=1,double parinterp=0,
                double morphedhistnorm=-1,
                int idebug=0)  const
  {

    //This is Alex Read's implementation of TH1D morphing.

    const char chname[50]="TH1D-interpolated";
    const char chtitle[50]="Interpolated histogram";

    // Return right away if one of the input histograms doesn't exist.
    if(!hist1) {
      cout << "ERROR! th1dmorph says first input histogram doesn't exist." << endl;
      return(0);
    }
    if(!hist2) {
      cout << "ERROR! th1dmorph says second input histogram doesn't exist." << endl;
      return(0);
    }

    // Extract bin parameters of input histograms 1 and 2. We haven't implemented
    // nonuniform binning (yet?).

    int nb1 = hist1->GetNbinsX();
    int nb2 = hist2->GetNbinsX();
    double xmin1 = hist1->GetXaxis()->GetXmin();
    double xmin2 = hist2->GetXaxis()->GetXmin();
    double xmax1 = hist1->GetXaxis()->GetXmax();
    double xmax2 = hist2->GetXaxis()->GetXmax();
    if (idebug > 0) {
      cout << nb1 << " " << xmin1 << " " << xmax1 << endl;
      cout << nb2 << " " << xmin2 << " " << xmax2 << endl;
    }

    // ......The weights (wt1,wt2) are the complements of the "distances" between
    //       the values of the parameters at the histograms and the desired
    //       interpolation point. For example, wt1=0, wt2=1 means that the
    //       interpolated histogram should be identical to input histogram 2.
    //       Check that they make sense. If par1=par2 then we can choose any
    //       valid set of wt1,wt2 so why not take the average?

    double wt1,wt2;
    if (par2 != par1) {
      wt1 = 1. - (parinterp-par1)/(par2-par1);
      wt2 = 1. + (parinterp-par2)/(par2-par1);
    }
    else {
      wt1 = 0.5;
      wt2 = 0.5;
    }

    //......Give a warning if this is an extrapolation.

    if (wt1 < 0 || wt1 > 1. || wt2 < 0. || wt2 > 1. || fabs(1-(wt1+wt2))
        > 1.0e-4) {
      cout << "Warning! th1dmorph: This is an extrapolation!! Weights are "
           << wt1 << " and " << wt2 << " (sum=" << wt1+wt2 << ")" << endl;
      //2014 - 08 - 18 Gianipez and P. Murat introduced the following changes for fixing an error
      // returned during the grid processing
      printf("parinterp = %5.5e par1 = %5.5e par2 = %5.5e \n", parinterp, par1, par2);
      if( wt1 <= 0.) {
        wt1 = 0.;
        wt2 = 1.;
      }else if(wt2 <= 0.) {
        wt2 = 0.;
        wt1 = 1.;
      }
      printf("Set weights to wt1 = %.5f wt2 = %.5f\n", wt1, wt2);
      //------------------------------------------------------------//
    }
    if (idebug >= 1) cout << "th1dmorph - Weights: " << wt1 << " " << wt2 << endl;

    //
    //......Perform the interpolation of histogram bin parameters. Use
    //      assignments instead of computation when input binnings
    //      are identical to assure best possible precision.

    double xminn=-1,xmaxn=-1;
    int nbn=0;
    double wtmin;

    wtmin = wt1; if (wt2 < wt1) wtmin = wt2;

    if (wtmin >= 0) {
      if (xmin1 == xmin2) {
        xminn = xmin1;
      } else {
        xminn = wt1*xmin1 + wt2*xmin2;
      }
      if (xmax1 == xmax2) {
        xmaxn = xmax1;
      } else {
        xmaxn = wt1*xmax1 + wt2*xmax2;
      }
      if (nb1 == nb2) {
        nbn = nb1;
      } else {
        nbn   = wt1*nb1   + wt2*nb2;
      }
    }
    //......If one of the weights is zero, then use the binnings of the
    //      histogram with nonzero weight.
    //      but reasonable with the histogram bin parameters.
    else {
      if (wt1 == 0) {
        xminn = xmin2; xmaxn = xmax2; nbn = nb2;
      } else if (wt2 == 0) {
        xminn = xmin1; xmaxn = xmax1; nbn = nb1;
      }
    }
    if (idebug >= 1) cout << "New hist: " << nbn << " " << xminn << " "
                          << xmaxn << endl;

    // Treatment for empty histograms: Return an empty histogram
    // with interpolated bins.

    if (hist1->GetSum() <= 0 || hist2->GetSum() <=0 ) {
      cout << "Warning! th1morph detects an empty input histogram. Empty interpolated histogram returned: "
           <<endl << "         " << chname << " - " << chtitle << endl;
      //      TH1D *morphedhist = (TH1D *)gROOT->FindObject(chname);
      //      if (morphedhist) delete morphedhist;
      TH1D* morphedhist = new TH1D(chname,chtitle,nbn,xminn,xmaxn);
      return(morphedhist);
    }
    if (idebug >= 1) cout << "Input histogram content sums: "
                          << hist1->GetSum() << " " << hist2->GetSum() << endl;
    // *
    // *......Extract the single precision histograms into double precision arrays
    // *      for the interpolation computation. The offset is because sigdis(i)
    // *      describes edge i (there are nbins+1 of them) while dist1/2
    // *      describe bin i. Be careful, ROOT does not use C++ convention to
    // *      number bins: dist1[ibin] is content of bin ibin where ibin runs from
    // *      1 to nbins. We allocate some extra space for the derived distributions
    // *      because there may be as many as nb1+nb2+2 edges in the intermediate
    // *      interpolated cdf described by xdisn[i] (position of edge i) and
    // *      sigdisn[i] (cummulative probability up this edge) before we project
    // *      into the final binning.

    double *dist1   = hist1->GetArray();
    double *dist2   = hist2->GetArray();
    double *sigdis1 = new double[1+nb1];
    double *sigdis2 = new double[1+nb2];
    double *sigdisn = new double[2+nb1+nb2];
    double *xdisn   = new double[2+nb1+nb2];
    double *sigdisf = new double[nbn+1];

    for(int i=0;i<2+nb1+nb2;i++) xdisn[i] = 0; // Start with empty edges
    sigdis1[0] = 0; sigdis2[0] = 0; // Start with cdf=0 at left edge

    for(int i=1;i<nb1+1;i++) {   // Remember, bin i has edges at i-1 and
      sigdis1[i] = dist1[i];       // i and i runs from 1 to nb.
    }
    for(int i=1;i<nb2+1;i++) {
      sigdis2[i] = dist2[i];
    }

    if (idebug >= 3) {
      for(int i=0;i<nb1+1;i++) {
        cout << i << " dist1" << dist1[i] << endl;
      }
      for(int i=0;i<nb1+1;i++) {
        cout << i << " dist2" << dist1[i] << endl;
      }
    }

    //......Normalize the distributions to 1 to obtain pdf's and integrate
    //      (sum) to obtain cdf's.

    double total = 0, norm1, norm2;
    for(int i=0;i<nb1+1;i++) {
      total += sigdis1[i];
    }
    if (idebug >=1) cout << "Total histogram 1: " <<  total << endl;
    for(int i=1;i<nb1+1;i++) {
      sigdis1[i] = sigdis1[i]/total + sigdis1[i-1];
    }
    norm1 = total;

    total = 0.;
    for(int i=0;i<nb2+1;i++) {
      total += sigdis2[i];
    }
    if (idebug >=1) cout << "Total histogram 22: " <<  total << endl;
    for(int i=1;i<nb2+1;i++) {
      sigdis2[i] = sigdis2[i]/total + sigdis2[i-1];
    }
    norm2 = total;

    // *
    // *......We are going to step through all the edges of both input
    // *      cdf's ordered by increasing value of y. We start at the
    // *      lower edge, but first we should identify the upper ends of the
    // *      curves. These (ixl1, ixl2) are the first point in each cdf from
    // *      above that has the same integral as the last edge.
    // *

    int ix1l = nb1;
    int ix2l = nb2;
    while(sigdis1[ix1l-1] >= sigdis1[ix1l]) {
      ix1l = ix1l - 1;
    }
    while(sigdis2[ix2l-1] >= sigdis2[ix2l]) {
      ix2l = ix2l - 1;
    }

    // *
    // *......Step up to the beginnings of the curves. These (ix1, ix2) are the
    // *      first non-zero points from below.

    int ix1 = -1;
    do {
      ix1 = ix1 + 1;
    } while(sigdis1[ix1+1] <= sigdis1[0]);

    int ix2 = -1;
    do {
      ix2 = ix2 + 1;
    } while(sigdis2[ix2+1] <= sigdis2[0]);

    if (idebug >= 1) {
      cout << "First and last edge of hist1: " << ix1 << " " << ix1l << endl;
      cout << "   " << sigdis1[ix1] << " " << sigdis1[ix1+1] << endl;
      cout << "First and last edge of hist2: " << ix2 << " " << ix2l << endl;
      cout << "   " << sigdis2[ix2] << " " << sigdis2[ix2+1] << endl;
    }

    //.......Need bin widths

    double dx1=(xmax1-xmin1)/double(nb1);
    double dx2=(xmax2-xmin2)/double(nb2);
    double dx=(xmaxn-xminn)/double(nbn);

    //......The first interpolated point should be computed now.

    int nx3 = 0;
    double x1,x2,x;
    x1 = xmin1 + double(ix1)*dx1;
    x2 = xmin2 + double(ix2)*dx2;
    x = wt1*x1 + wt2*x2;
    xdisn[nx3] = x;
    sigdisn[nx3] = 0;
    if(idebug >= 1) {
      cout << "First interpolated point: " << xdisn[nx3] << " "
           << sigdisn[nx3] << endl;
      cout << "                          " << x1 << " <= " << x << " <= "
           << x2 << endl;
    }

    //......Loop over the remaining point in both curves. Getting the last
    //      points may be a bit tricky due to limited floating point
    //      precision.

    if (idebug >= 1) {
      cout << "----BEFORE while with ix1=" << ix1 << ", ix1l=" << ix1l
           << ", ix2=" << ix2 << ", ix2l=" << ix2l << endl;
      for(int i=ix1;i<=ix1l;i++) {cout << "   1: " << i << " " << sigdis1[i] << endl;}
      for(int i=ix2;i<=ix2l;i++) {cout << "   2: " << i << " " << sigdis2[i] << endl;}
    }

    double yprev = -1; // The probability y of the previous point, it will
    //get updated and used in the loop.
    double y,x20,x21,y20,y21; // Interpolation points along cdfs 0,1,2
    double x10,x11,y10,y11;

    while((ix1 < ix1l) || (ix2 < ix2l)) {
      if (idebug >= 1 ) cout << "----Top of while with ix1=" << ix1
                             << ", ix1l=" << ix1l << ", ix2=" << ix2
                             << ", ix2l=" << ix2l << endl;

      //......Increment to the next lowest point. Step up to the next
      //      kink in case there are several empty (flat in the integral)
      //      bins.

      int i12type = -1; // Tells which input distribution we need to
      // see next point of.

      if (((sigdis1[ix1+1] <= sigdis2[ix2+1]) || (ix2 == ix2l)) && ix1 < ix1l) {
        ix1 = ix1 + 1;
        // try to fix empty bin holes!!!!
        //      while(sigdis1[ix1+1] <= sigdis1[ix1] && ix1 < ix1l) {
        //         ix1 = ix1 + 1;
        //      }
        //empty bin fix??? while(sigdis1[ix1+1] <= sigdis1[ix1] && ix1 < ix1l) {
        while(sigdis1[ix1+1] < sigdis1[ix1] && ix1 < ix1l) {
          ix1 = ix1 + 1;
        }
        i12type = 1;
      } else if (ix2 < ix2l) {
        ix2 = ix2 + 1;
        //empty bin fix ?? while(sigdis2[ix2+1] <= sigdis2[ix2] && ix2 < ix2l) {
        while(sigdis2[ix2+1] < sigdis2[ix2] && ix2 < ix2l) {
          ix2 = ix2 + 1;
        }
        i12type = 2;
      }
      if (i12type == 1) {
        if (idebug >= 3) {
          cout << "Pair for i12type=1: " << ix1 << " " << ix2 << " " << sigdis2[ix2] << " "
               << sigdis1[ix1] << " " << sigdis2[ix2+1] << endl;
        }
        x1  = xmin1 + double(ix1)*dx1 ;
        y   = sigdis1[ix1];
        x20 = double(ix2)*dx2 + xmin2;
        x21 = x20 + dx2;
        y20 = sigdis2[ix2];
        y21 = sigdis2[ix2+1];

        //......Calculate where the cummulative probability y in distribution 1
        //      intersects between the 2 points from distribution 2 which
        //      bracket it.

        if (y21 > y20) {
          x2 = x20 + (x21-x20)*(y-y20)/(y21-y20);
        }
        else {
          x2 = x20;
        }
      } else {
        if (idebug >= 3) {
          cout << "Pair for i12type=2: " << sigdis1[ix1] << " " << sigdis2[ix2]
               << " " << sigdis1[ix1+1] << endl;
        }
        x2 = xmin2 + double(ix2)*dx2 ;
        y = sigdis2[ix2];
        x10 = double(ix1)*dx1 + xmin1;
        x11 = x10 + dx1;
        y10 = sigdis1[ix1];
        y11 = sigdis1[ix1+1];

        //......Calculate where the cummulative probability y in distribution 2
        //      intersects between the 2 points from distribution 1 which
        //      brackets it.

        if (y11 > y10) {
          x1 = x10 + (x11-x10)*(y-y10)/(y11-y10);
        } else {
          x1 = x10;
        }
      }

      //......Interpolate between the x's in the 2 distributions at the
      //      cummulative probability y. Store the (x,y) for provisional
      //      edge nx3 in (xdisn[nx3],sigdisn[nx3]). nx3 grows for each point
      //      we add the the arrays. Note: Should probably turn the pair into
      //      a structure to make the code more object-oriented and readable.

      x = wt1*x1 + wt2*x2;
      if (y >= yprev) { // bugfix for empty bins?!?!?!?!
        nx3 = nx3+1;
        if (idebug >= 1) {
          cout << " ---> y > yprev: i12type=" << i12type << ", nx3="
               << nx3 << ", x= " << x << ", y=" << y << ", yprev=" << yprev
               << endl;
        }
        yprev = y;
        xdisn[nx3] = x;
        sigdisn[nx3] = y;
        if(idebug >= 1) {
          cout << "    ix1=" << ix1 << ", ix2= " << ix2 << ", i12type= "
               << i12type << ", sigdis1[ix1]=" << sigdis1[ix1] << endl;
          cout << "        " << ", nx3=" << nx3 << ", x=" << x << ", y= "
               << sigdisn[nx3] << endl;
        }
      }
    }
    if (idebug >=3) for (int i=0;i<=nx3;i++) {
        cout << " nx " << i << " " << xdisn[i] << " " << sigdisn[i] << endl;
      }

    // *......Now we loop over the edges of the bins of the interpolated
    // *      histogram and find out where the interpolated cdf 3
    // *      crosses them. This projection defines the result and will
    // *      be stored (after differention and renormalization) in the
    // *      output histogram.
    // *
    // *......We set all the bins following the final edge to the value
    // *      of the final edge.

    x = xminn + double(nbn)*dx;
    int ix = nbn;

    if (idebug >= 1) cout << "------> Any final bins to set? " << x << " "
                          << xdisn[nx3] << endl;
    while(x >= xdisn[nx3]) {
      sigdisf[ix] = sigdisn[nx3];
      if (idebug >= 2) cout << "   Setting final bins" << ix << " " << x
                            << " " << sigdisf[ix] << endl;
      ix = ix-1;
      x = xminn + double(ix)*dx;
    }
    int ixl = ix + 1;
    if (idebug >= 1) cout << " Now ixl=" << ixl << " ix=" << ix << endl;

    // *
    // *......The beginning may be empty, so we have to step up to the first
    // *      edge where the result is nonzero. We zero the bins which have
    // *      and upper (!) edge which is below the first point of the
    // *      cummulative distribution we are going to project to this
    // *      output histogram binning.
    // *

    ix = 0;
    x = xminn + double(ix+1)*dx;
    if (idebug >= 1) cout << "Start setting initial bins at x=" << x << endl;
    while(x <= xdisn[0]) {
      sigdisf[ix] = sigdisn[0];
      if (idebug >= 1) cout << "   Setting initial bins " << ix << " " << x
                            << " " << xdisn[1] << " " << sigdisf[ix] << endl;
      ix = ix+1;
      x = xminn + double(ix+1)*dx;
    }
    int ixf = ix;

    if (idebug >= 1)
      cout << "Bins left to loop over:" << ixf << "-" << ixl << endl;

    // *......Also the end (from y to 1.0) often comes before the last edge
    // *      so we have to set the following to 1.0 as well.

    int ix3 = 0; // Problems with initial edge!!!
    for(ix=ixf;ix<ixl;ix++) {
      x = xminn + double(ix)*dx;
      if (x < xdisn[0]) {
        y = 0;
      } else if (x > xdisn[nx3]) {
        y = 1.;
      } else {
        while(xdisn[ix3+1] <= x && ix3 < 2*nbn) {
          ix3 = ix3 + 1;
        }
        if (xdisn[ix3+1]-x > 1.1*dx2) { // Empty bin treatment
          y = sigdisn[ix3+1]; //y = sigdisn[ix3+1];
          if(idebug>=1) cout << "Empty bin treatment " << ix3+1 << " " << y << " " << sigdisn[ix3+2] << " " << sigdisn[ix3-1] << endl;
        }
        else if (xdisn[ix3+1] > xdisn[ix3]) { // Normal bins
          y = sigdisn[ix3] + (sigdisn[ix3+1]-sigdisn[ix3])
            *(x-xdisn[ix3])/(xdisn[ix3+1]-xdisn[ix3]);
        } else {  // Is this ever used?
          y = 0;
          cout << "Warning - th1dmorph: This probably shoudn't happen! "
               << endl;
          cout << "Warning - th1dmorph: Zero slope solving x(y)" << endl;
        }
      }
      sigdisf[ix] = y;
      if (idebug >= 3) {
        cout << ix << ", ix3=" << ix3 << ", xdisn=" << xdisn[ix3] << ", x="
             << x << ", next xdisn=" << xdisn[ix3+1] << endl;
        cout << "   cdf n=" << sigdisn[ix3] << ", y=" << y << ", next point="
             << sigdisn[ix3+1] << endl;
      }
    }

    //......Differentiate interpolated cdf and return renormalized result in
    //      new histogram.

    //    TH1D *morphedhist = (TH1D *)gROOT->FindObject(chname);
    //    if (morphedhist) delete morphedhist;
    TH1D* morphedhist = new TH1D(chname,chtitle,nbn,xminn,xmaxn);

    double norm = morphedhistnorm;
    // norm1, norm2, wt1, wt2 are computed before the interpolation
    if (norm <= 0) {
      if (norm1 == norm2) {
        norm = norm1;
      } else {
        norm   = wt1*norm1   + wt2*norm2;
      }
    }

    for(int ixx=nbn-1;ixx>-1;ixx--) {
      x = xminn + double(ixx)*dx;
      y =  sigdisf[ixx+1]-sigdisf[ixx];
      if (y<0) cout << "huh??? " << ixx << " " << sigdisf[ixx] << " " << sigdisf[ixx+1] << endl;
      morphedhist->SetBinContent(ixx+1,y*norm);
    }

    //......Clean up the temporary arrays we allocated.

    delete [] sigdis1; delete [] sigdis2;
    delete [] sigdisn; delete [] xdisn; delete [] sigdisf;

    //......All done, return the result.
    //    morphedhist->Draw("same");
    return(morphedhist);
  }

}
