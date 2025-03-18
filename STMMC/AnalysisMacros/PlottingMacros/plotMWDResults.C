// Generates the analysis path of the MWD algorithm
// See plotMWDResults.sh for usage example
// Original author - Pawel Plesniak

#include <limits>

void customErrorHandler(int level, Bool_t abort, const char* location, const char* message) {
    /*
        Description
            Define a custom error handler that won't print the stack trace but will print an error message and exit.
    */
    std::cerr << message << std::endl;
    if (level > kInfo)
        exit(1);
};

void collectData(const std::string fileName, const std::string treeName, std::vector<double> &times, std::vector<double> &energies) {
    /*
        Description
            Collects all the required data from virtual detector TTrees

        Arguments
            fileName - as documented in function "plotMWDResults"
            treeName - as documented in function "plotMWDResults"
            times - as documented in function "plotMWDResults"
            energies - as documented in function "plotMWDResults"

        Variables
            file - ROOT TFile interface
            branches - ROOT TFile interface to TTree branches
            branchNames - vector of branch names
            dataTime - time from file
            dataE - energy from file
            entries - number of entries in the TTree
    */
    std::cout << "Processing file " << fileName << std::endl;

    // Get the branch
    TFile *file = new TFile(fileName.c_str());
    if (!file || file->IsZombie()) {
        Fatal("collectData", "Failed to open the file.");
    };
    TTree *tree = (TTree*)file->Get(treeName.c_str());
    if (!tree)
        Fatal("collectData", "Requested tree does not exist in the file.");

    // Get the list of branches to check if they exist
    TObjArray *branches = tree->GetListOfBranches();
    std::vector<std::string> branchNames;
    for (int i = 0; i < branches->GetEntries(); ++i) {
        TBranch *branch = dynamic_cast<TBranch*>(branches->At(i));
        if (branch)
            branchNames.push_back(branch->GetName());
    };

    // If the branches exist, assign them to the appropriate variables
    double dataE;
    uint32_t dataTime;
    if (std::find(branchNames.begin(), branchNames.end(), "time") != branchNames.end())
        tree->SetBranchAddress("time", &dataTime);
    else
        Fatal("collectData", "Requested branch 'time' does not exist in the file.");
    if (std::find(branchNames.begin(), branchNames.end(), "E") != branchNames.end())
        tree->SetBranchAddress("E", &dataE);
    else
        Fatal("collectData", "Requested branch 'E' does not exist in the file.");

    // Get the number of entries
    int entries = tree->GetEntries();

    // Set up constant for converting ADC time from ADC clock ticks to real time
    const double tADC = 3.125;

    // Collect the data
    for (int i = 0; i < entries; i++) {
        // Get the corresponding entry
        tree->GetEntry(i);

        // Collect the data
        times.push_back(dataTime * tADC);
        energies.push_back(dataE);
    };
    // Check if data has been collected
    if (times.size() == 0 || energies.size() == 0)
        Fatal("collectData", "No data was collected from this file");
    if (times.size() != energies.size())
        Fatal("collectData", "Unequal number of data points were collected");

    // Clear
    file->Close();
    delete file;
    std::cout << "Finished processing file " << fileName << "\n" << std::endl;
    return;
};


std::string convertBinWidthToStr(double binWidth) {
    /*
        Description
            Converts the bin width as a double to a string and returns it

        Arguments
            binWidth - bin width

        Variables
            stream - input string stream to read the bin width
            binWidthStr - bin width as a str
    */
    // Set up the stream to read in the bin width
    std::stringstream stream;

    // Set the precision and read it in. If the bin width is an integer, set 0 decimal places, otherwise set 3
    stream << std::fixed << std::setprecision(((binWidth - (int)binWidth) < std::numeric_limits<double>::epsilon()) ? 0 : 3) << binWidth;

    // Convert the stream object to a string
    std::string binWidthStr = stream.str();

    // Clear
    stream.str("");
    std::cout << std::defaultfloat;
    return binWidthStr;
};

void plot(std::vector<double> times, std::vector<double> energies, std::vector<double> binWidths, double ERed, bool convertkeVToMeV, const double signalAcceptance, bool highResolution, bool timeCuts, bool flashCut347) {
    /*
        Description
            Plots the full, reduced, and signal spectra

        Arguments
            times - as documented in function "plotMWDResults"
            energies - as documented in function "plotMWDResults"
            binWidths - as documented in function "plotMWDResults"
            ERed - as documented in function "plotMWDResults"
            convertkeVToMeV - as documented in function "plotMWDResults"
            signalAcceptance - as documented in function "plotMWDResults"
            highResolution - as documented in function "plotMWDResults"
            timeCuts - as documented in function "plotMWDResults"
            flashCut347 - as documented in function "plotMWDResults"

        Variables
            eMax - maximum energy in plotting data
            eMin - minomum energy in plotting data
            eRangeFull - energy range of data to be plotted
            eRangeRed - energy range of the reduced spectrum
            e347 - energy of 347 keV signal
            e844 - energy of 844 keV signal
            e1809 - energy of 1809 keV signal
            eRange347 - energy range of 347 keV signal
            eRange844 - energy range of 844 keV signal
            eRange1809 - energy range of 1809 keV signal
            eMin347 - minimum of 347 keV signal plot
            eMin844 - minimum of 844 keV signal plot
            eMin1809 - minimum of 1809 keV signal plot
            eMax347 - maximum of 347 keV signal plot
            eMax844 - maximum of 844 keV signal plot
            eMax1809 - maximum of 1809 keV signal plot
            tMicrospill - microspill duration in ns
            tMin347 - minimum time for 347keV signal acceptance
            tMax347 - maximum time for 347keV signal acceptance
            tMod347 - time modulus for 347keV signal acceptance
            tMin844 - minimum time for 844keV signal acceptance
            tMax844 - maximum time for 844keV signal acceptance
            tMod844 - time modulus for 844keV signal acceptance
            tMin1809 - minimum time for 1809keV signal acceptance
            tMax1809 - maximum time for 1809keV signal acceptance
            tMod1809 - time modulus for 1809keV signal acceptance
            vars - temporary buffer to store data if it requires unit conversion from keV to MeV
            order - order of plot generation histograms
            eMinMax - energy ranges of histograms as {minimum, maximum}
            eRanges - energy ranges of the histograms as {maximum - minimum}
            nOrder - number of plot categories
            nBins - bin numbers in histograms
            binWidthStrs - bin widths as strings for histograms
            unit - energy unit as either keV or MeV
            sigRange - controls when the signal specific plot title extensions are added
            fileNames - vector of plot file names
            titles - vector of plot titles
            px - number of x pixels for canvases
            py - number of y pixels for canvases
            canvases - vector of all canvases
            plot - order iterator
            lowResAxisTextSize - low resolution axis text size
            c - canvases iterator
            hists - vector of histograms
            nEntries - number of data points to use with histograms
            e - energy buffer vector
            t - time buffer vector
            lx1 - defines an x coordinate for stat boxes and legend
            ly1 - defines a y coordinate for stat boxes and legend
            lx2 - defines the other x coordinate for stat boxes and legend
            ly2 - defines the other y coordinate for stat boxes and legend
            count - buffer vector for histogram bin checks

    */
    std::cout << std::string("Generating plots with ")  + (highResolution ? "high resolution, " : "low resolution, ") + (timeCuts ? "time cuts, " : "no time cuts, ") + (convertkeVToMeV ? "in MeV" : "in MeV") << std::endl;

    // Define the energy parameters in keV
    double eMax = *std::max_element(energies.begin(), energies.end()) + 1, eMin = *std::min_element(energies.begin(), energies.end());
    if (eMin < 0) {
        std::cout << "==========Warning==========" << std::endl;
        std::cout << "Negative minimum energy of " << *std::min_element(energies.begin(), energies.end()) << " keV" << std::endl;
        std::cout << "===========================" << std::endl;
    };
    double eRangeFull = eMax - eMin;
    double eRangeRed  = ERed - eMin;
    double e347  = 347, eRange347  = e347  * signalAcceptance,  eMin347  = e347  - eRange347,   eMax347  = e347  + eRange347;
    double e844  = 844, eRange844  = e844  * signalAcceptance,  eMin844  = e844  - eRange844,   eMax844  = e844  + eRange844;
    double e1809 = 809, eRange1809 = e1809 * signalAcceptance,  eMin1809 = e1809 - eRange1809,  eMax1809 = e1809 + eRange1809;

    // Double the energy range to include both sides of the acceptance
    eRange347  *=2;
    eRange844  *=2;
    eRange1809 *=2;

    // Define the time parameters in ns
    const double tMicroSpill = 1695;
    const double tMin347  = flashCut347 ? 200 : 300,    tMax347  = 700,     tMod347  = tMicroSpill;
    const double tMin844  = 492000,                     tMax844  = 1330000, tMod844  = tMax844;
    const double tMin1809 = 500,                        tMax1809 = 1600,    tMod1809 = tMicroSpill;

    // Convert everything to keV if needed
    double* vars[] = {&e347,  &eMin347,  &eMax347,  &eRange347,
                      &e844,  &eMin844,  &eMax844,  &eRange844,
                      &e1809, &eMin1809, &eMax1809, &eRange1809,
                      &eMin,  &eMax,     &ERed,     &eRangeFull, &eRangeRed};
    if (convertkeVToMeV) {
        double conversion = 1000;
        std::transform(std::begin(vars),  std::end(vars),  std::begin(vars),  [](double* x)-> double*{ *x /= 1000; return x;});
        std::transform(energies.begin(),  energies.end(),  energies.begin(),  [conversion](double x) { return x / conversion; });
        std::transform(binWidths.begin(), binWidths.end(), binWidths.begin(), [conversion](double x) { return x / conversion; });
    };
    std::vector<std::string>            order   = {"Full",          "Reduced",      "347",              "844",              "1809"};
    std::vector<std::vector<double>>    eMinMax = {{eMin, eMax},    {eMin, ERed},   {eMin347, eMax347}, {eMin844, eMax844}, {eMin1809, eMax1809}};  // Energy ranges for TH1Ds
    std::vector<double>                 eRanges = {eRangeFull,      eRangeRed,      eRange347,          eRange844,          eRange1809};
    const int nOrder = order.size();
    std::vector<int> nBins(nOrder);
    std::transform(eRanges.begin(), eRanges.end(), binWidths.begin(), nBins.begin(), [](double E, double BW) {if (BW == 0) throw std::runtime_error("Illegal zero bin width, exiting."); return E / BW;});

    // Convert non str types to str types
    std::vector<std::string> binWidthStrs;
    for (double binWidth : binWidths)
        binWidthStrs.push_back(convertBinWidthToStr(binWidth));
    std::string unit = convertkeVToMeV ? "MeV" : "keV";

    // Construct file names
    const int sigRange = 1;
    std::vector<std::string> fileNames;
    std::string signalExtension = std::string("") + (timeCuts ? ".time-cut" : ".no-time-cut") + (flashCut347 ? ".347-flash" : "");
    for (int i = 0; i < nOrder; i++)
        fileNames.push_back("MWD." + order[i] + "." + unit + "." + (highResolution ? "high" : "low") + (i > sigRange ? signalExtension : "") + ".png");

    // Construct plot titles
    std::vector<std::string> titles;
    for (int i = 0; i < nOrder; i++)
        titles.push_back("Reconstructed signal energy; Energy [" + unit + "]; Count / " + binWidthStrs[i] + " " + unit);

    // Set up the TCanvases
    const int px = (highResolution ? 1500 : 750), py = (highResolution ? 1000 : 500);
    std::vector<TCanvas*> canvases;
    for (std::string plot : order)
        canvases.push_back(new TCanvas(("c" + plot).c_str(), ("c" + plot).c_str(), px, py));

    // Apply TCanvas formatting
    const double lowResAxisTextSize = 0.06;
    gStyle->SetTitleFontSize(lowResAxisTextSize);
    gStyle->SetTitleAlign(33);
    gStyle->SetTitleX(0.75);
    if (!highResolution) {
        gStyle->SetTitleX(0.99);
        for (TCanvas* c : canvases) {
            c->SetBottomMargin(0.14);
            c->SetLeftMargin(0.175);
            c->Update();
        };
    };

    // Set up TH1Ds
    std::vector<TH1D*> hists;
    for (int i = 0; i < nOrder; i++)
        hists.emplace_back(new TH1D(("h" + order[i]).c_str(), ("h" + order[i]).c_str(), nBins[i], eMinMax[i][0], eMinMax[i][1]));

    // Populate the TH1Ds
    const int nEntries = energies.size();
    double e = 0.0, t = 0.0;
    for (uint i = 0; i < nEntries; i++) {
        e = energies[i];
        t = times[i];
        hists[0]->Fill(e);
        if (e < ERed)
            hists[1]->Fill(e);
        if ((!timeCuts && eMin347 < e  && e < eMax347)  || (timeCuts && tMin347 < fmod(t, tMod347)   && fmod(t, tMod347) < tMax347   && eMin347 < e  && e < eMax347))
            hists[2]->Fill(e);
        if ((!timeCuts && eMin844 < e  && e < eMax844)  || (timeCuts && tMin844 < fmod(t, tMod844)   && fmod(t, tMod844) < tMax844   && eMin844 < e  && e < eMax844))
            hists[3]->Fill(e);
        if ((!timeCuts && eMin1809 < e && e < eMax1809) || (timeCuts && tMin1809 < fmod(t, tMod1809) && fmod(t, tMod1809) < tMax1809 && eMin1809 < e && e < eMax1809))
            hists[4]->Fill(e);
    };

    // Set up the stat boxes
    const double lx1 = (highResolution ? 0.7 : 0.6);
    const double ly1 = lx1;
    const double lx2 = 0.9;
    const double ly2 = 0.9;
    gStyle->SetStatX(0.9);
    gStyle->SetStatY(0.9);
    gStyle->SetStatW(lx2 - lx1);
    gStyle->SetStatH(ly2 - ly1);

    // Check the hists for overflow and underflow bins
    int count = 0;
    for (int i = 0; i < nOrder; i++) {
        count = hists[i]->GetBinContent(0);
        if (count) {
            std::cout << "In plot " << fileNames[i] << ", number of underflow bins is " << count << std::endl;
            std::cout << "This histogram's lower range is: " << hists[i]->GetXaxis()->GetXmin() << std::endl;
            Fatal("makePlot", "Generated plot has values in the underflow bins");
        };
        count = hists[i]->GetBinContent(hists[i]->GetNbinsX() + 1);
        if (count) {
            std::cout << "In plot " << fileNames[i] << ", number of overflow bins is " << count << std::endl;
            std::cout << "This histogram's upper range is: " << hists[i]->GetXaxis()->GetXmax() << std::endl;
            Fatal("makePlot", "Generated plot has values in the overflow bins");
        };
        count = hists[i]->GetEntries() - hists[i]->Integral(1, hists[i]->GetNbinsX());
        if (count) { // scaling affects integral but not entries
            std::cout << "In plot " << fileNames[i] << ", difference between entries and integral is " << count << std::endl;
            Fatal("makePlot", "Generated plot has values in the underflow/overflow bins");
        };
    };

    // Draw and save the unstacked plots
    for (int i = 0; i < nOrder; i++) {
        canvases[i]->cd();
        count = hists[i]->GetEntries();
        if (count == 0) {
            std::cout << fileNames[i] << " has no entries, not saving" << std::endl;
            continue;
        };
        hists[i]->SetTitle(titles[i].c_str());
        hists[i]->SetMinimum(0); // Sets the y minimum
        hists[i]->Draw("HIST");
        if (!highResolution) {
            hists[i]->GetXaxis()->SetLabelSize(lowResAxisTextSize);
            hists[i]->GetXaxis()->SetTitleSize(lowResAxisTextSize);
            hists[i]->GetYaxis()->SetLabelSize(lowResAxisTextSize);
            hists[i]->GetYaxis()->SetTitleSize(lowResAxisTextSize);
        };
        canvases[i]->SaveAs(fileNames[i].c_str());
    };

    // Clean up from the unstacked plots
    for (int i = nOrder - 1; i > -1; i--) {
        hists[i]->Delete();
        canvases[i]->Close();
    };

    // Done
    std::cout << "Finished\n" << std::endl;
    return;
};

void makePlots(std::vector<double> &times, std::vector<double> &energies, double ERed, const double signalAcceptance, std::vector<double> &binWidths) {
    /*
        Description
            Sets up loops for function "plot"

        Arguments
            times - as documented in function "plotMWDResults"
            energies - as documented in function "plotMWDResults"
            ERed - as documented in function "plotMWDResults"
            binWidths - as documented in function "plotMWDResults"

        Variables
    */

    std::vector<bool> bools = {true, false};
    for (bool highResolution : bools) {
        for (bool convertkeVToMeV : bools) {
            plot(times, energies, binWidths, ERed, convertkeVToMeV, signalAcceptance, highResolution, false, false);
            plot(times, energies, binWidths, ERed, convertkeVToMeV, signalAcceptance, highResolution, true,  false);
            plot(times, energies, binWidths, ERed, convertkeVToMeV, signalAcceptance, highResolution, true,  true);
        };
    };

    return;
};


void plotMWDResults(std::string fileName, std::string treeName, double ERed = 2000.0, const double signalAcceptance = 0.1, double binWidthFull = 500, double binWidthRed = 10, double binWidth347 = 2, double binWidth844 = 2, double binWidth1809 = 2) {
    /*
        Description
            Plots spectra measured by the detectors using the MWD algorithm. Names the files as
                MWD.<full/red/signal>.<keV/MeV>.<high/low>.<time-cut/no-time-cut>.<347-flash/>.png
            such that
                <full/red/signal> - determines the range of the plot, with "full" plotting the full spectrum, "red" plotting the spectrum up to ERed, and "signal" corresponding to each signal in the STM energy range
                <keV/MeV> - units of the energy (x) axis
                <high/low> - whether the plot is in high or low resolution
                <time-cut/no-time-cut> - defines whether or not the relevant time cut has been applied
                <347-flash/> - if present, uses the 347keV flash cut instead of the off-flash cut

        Arguments
            fileName - name of file containing MWD derived spectra in a ROOT file as a relative path to cwd
            treeName - name of tree containing all the data

        Optional arguments
            ERed - reduced spectrum max energy [MeV]
            signalAcceptance - for the signal plots, this controls the width of the plotted region as a multiple of the signal energy
            binWidthFull - full energy spectrum bin width
            binWidthRed - reduced energy spectrum bin width
            binWidth347 - bin width of the 347keV signal plot
            binWidth844 - bin width of the 844keV signal plot
            binWidth1809 - bin width of the 1809keV signal plot

        Variables
            times - vector of times converted from ADC clock tick times to real times from fileName
            energies - vector of energies measured by detector [keV]
            binWidths - vector of bin widths
    */

    // Update global parameters
    SetErrorHandler(customErrorHandler);
    gROOT->SetBatch(kTRUE);

    // Set up the data storage variables
    std::vector<double> times, energies;

    // Collect the data from the file
    collectData(fileName, treeName, times, energies);

    // Gather the bin widths together
    std::vector<double> binWidths = {binWidthFull, binWidthRed, binWidth347, binWidth844, binWidth1809};

    // Generate the plots
    makePlots(times, energies, ERed, signalAcceptance, binWidths);

    return;
};
