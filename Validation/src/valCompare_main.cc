//
// Compare two files of histograms and report if 
// they are consistent.  For Validation procedures.
//

#include "TFile.h"

#include <iostream>
#include <string>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "Validation/inc/TValCompare.hh"

void valCompare_usage() {

  std::cout << 
"  \n"
"      valCompare [OPTIONS] FILE1 FILE2\n"
"  \n"
"  Compare two root histogram files.  Only the 1D histograms, \n"
"  profiles and efficiencies will be compared - other objects are skipped. \n"
"  The first histogram file on the command line appears as a histogram \n"
"  in the plots and first in listing or characteristics. \n"
"  Each comparison result is given a status:\n"
"  0 = perfect match\n"
"  1 = matches > 99.9% in K-S or fractional test\n"
"  2 = matches > 99% in K-S or fractional test\n"
"  10 = one or both histograms is empty\n"
"  11 = can't compare\n"
"  \n"
"  Examples:\n"
"    - browse plots which fail lose and tight criteria\n"
"    valCompare -b -l 3 FILE1 FILE2\n"
"    - print a grand summary of the comparison\n"
"    valCompare -s FILE1 FILE2\n"
"    - print a line for each plot failing exact match\n"
"    valCompare -r -l 1 FILE1 FILE2\n"
"    - put all comparisons into a pdf file with 2x2 on a page\n"
"    valCompare -2 -p result.pdf FILE1 FILE2\n"
"  \n"
"  -h print help\n"
"  -v INT verbose level (default=1)\n"
"  -l INT  select plots to show - lower limit to status (0-11)\n"
"  -g INT  select plots to show - upper limit to status (0-11)\n"
"  -a FLOAT  threshold for tight agreement (default = 0.999)\n"
"  -e FLOAT  threshold for loose agreement (default = 0.99)\n"
"  -i treat samples as statistcally independent instead of ~identical\n"
"  -c FLOAT  scaling for 1st file\n"
"  -d FLOAT  scaling for 2nd file\n"
"  -m INT mode: 0= no scaling, 1 = scale 2nd to 1st, 2=scale \n"
"      to input values for mode 2, switches and b are required. \n" 
"      Only effects fraction comparison.\n"
"  -b browse the plots\n"
"  -s print a summary\n"
"  -r print a one-line report for each histogram\n"
"  -1 put plots on page 1x2\n"
"  -2 put plots on page 2x2\n"
"  -u ignore underflows in comparison\n"
"  -o ignore overflows in comparison\n"
"  -p FILE  PDF file output like dir/results.pdf\n"
"  -w FILE  web page output like dir/dir/results.html\n"
	 << std::endl;
  return;
}



int main (int argc, char **argv)
{

  int verbose = 1;
  float loose = 9999.0;
  float tight = 9999.0;
  float scale1 = 1.0;
  float scale2 = 1.0;
  int mode = 0;
  bool qBrowse = false;
  bool qSummary = false;
  bool qReport = false;
  bool q12 = false;
  bool q22 = false;
  int llim = -1;
  int ulim = 999;
  int under = 1;
  int over  = 1;
  int indep = 0;

  char* webPage = nullptr;
  char* pdfFile = nullptr;

  char c;

  opterr = 0;
  while ((c = getopt (argc, argv, "hv:a:e:ibc:d:m:qsr12l:g:uo:p:w:")) != -1)
    switch (c)
      {
      case 'h':
        valCompare_usage();
        return 0;
      case 'v':
        verbose = atoi(optarg);
        break;

      case 'l':
        llim = atoi(optarg);
        break;
      case 'g':
        ulim = atoi(optarg);
        break;
      case 'a':
        llim = atoi(optarg);
        break;
      case 'e':
        ulim = atoi(optarg);
        break;
      case 'i':
        indep = 1;
        break;
      case 'c':
        scale1 = atof(optarg);
        break;
      case 'd':
        scale2 = atof(optarg);
        break;
      case 'm':
        mode = atoi(optarg);
        break;

      case 'b':
        qBrowse = true;
        break;
      case 's':
        qSummary = true;
        break;
      case 'r':
        qReport = true;
        break;
      case '1':
        q12 = true;
        break;
      case '2':
        q22 = true;
        break;
      case 'u':
        under = 1;
        break;
      case 'o':
        over = 1;
        break;
      case 'p':
        pdfFile = optarg;
        break;
      case 'w':
        webPage = optarg;
        break;
      case '?':
        valCompare_usage();
        return 1;
      default:
        valCompare_usage();
        return 1;
      }

  if( optind >= argc ) {
    printf("ERROR - need two histogram file names on the command line\n");
    valCompare_usage();
    return 1;
  }
  TValCompare pp;
  pp.SetVerbose(verbose);
  pp.SetFile1(argv[optind]);
  pp.SetFile2(argv[optind+1]);
  pp.SetMinStat(llim);
  pp.SetMaxStat(ulim);

  TValPar& pr = pp.GetPar();
  // once these are set, then indep won't override them
  // only set them if actually requested
  if(loose<9998.0) pr.SetLoose(loose);
  if(tight<9998.0) pr.SetTight(tight);
  pr.SetIndependent(indep);
  pr.SetUnder(under);
  pr.SetOver(over);
  pr.SetScale1(scale1);
  pr.SetScale2(scale2);
  pr.SetMode(mode);

  pp.Analyze();
  if(qReport) pp.Report();
  if(qSummary) pp.Summary();
  std::string opt;
  if(q12) {
    opt.append("1X2");
  } else if (q22) {
    opt.append("2X2");
  }
  if(qBrowse) pp.Display(opt.c_str());
  if(pdfFile) pp.SaveAs(pdfFile,opt.c_str());
  if(webPage) pp.SaveAs(webPage,opt.c_str());

  return 0;
}

