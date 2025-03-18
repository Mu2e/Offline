// Generates spectra of data generated from STM MC method, see doc db 51487 for more information
// See plot.sh for usage examples
// Note - it is recommended to run this on mu2ebuild02 as it takes a LOT of memory to run this macro
// Original author: Pawel Plesniak

#include <TError.h>
#include <iostream>
#include <limits>
#include <math.h>

void customErrorHandler(int level, Bool_t abort, const char* location, const char* message) {
    /*
        Description
            Define a custom error handler that won't print the stack trace but will print an error message and exit.
    */
    std::cerr << message << std::endl;
    if (level > kInfo)
        exit(1);
};

void collectVirtualdetectorData(const std::string fileName, const std::string treeName, std::vector<double> &energies, std::vector<double> &times, std::vector<ULong64_t> &virtualdetectorIds, std::vector<int> &pdgIds) {
    /*
        Description
            Collects all the required data from virtual detector TTrees

        Arguments
            fileName - as documented in function "plot"
            treeName - as documented in function "plot"
            energies - vector of energies associated with the relevant dataset, in MeV
            times - vector of times associated with the relevant dataset, in ns
            virtualdetectorIds - vector of virtual detector IDs associated with the relevant dataset
            pdgIds - vector of PDG IDs associated with the relevant dataset

        Variables
            file - ROOT TFile interface
            branches - ROOT TFile interface to TTree branches
            branchNames - vector of branch names
            dataVirtualDetectorId - virtual detector ID from file
            dataPdgId - PDG ID from file
            dataE - energy from file
            dataTime - time from file
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
    ULong64_t dataVirtualdetectorId;
    int dataPdgId;
    double dataE, dataTime;
    if (std::find(branchNames.begin(), branchNames.end(), "virtualdetectorId") != branchNames.end())
        tree->SetBranchAddress("virtualdetectorId", &dataVirtualdetectorId);
    else
        Fatal("collectData", "Requested branch 'virtualdetectorId' does not exist in the file.");
    if (std::find(branchNames.begin(), branchNames.end(), "pdgId") != branchNames.end())
        tree->SetBranchAddress("pdgId", &dataPdgId);
    else
        Fatal("collectData", "Requested branch 'pdgId' does not exist in the file.");
    if (std::find(branchNames.begin(), branchNames.end(), "E") != branchNames.end())
        tree->SetBranchAddress("E", &dataE);
    else
        Fatal("collectData", "Requested branch 'E' does not exist in the file.");
    if (std::find(branchNames.begin(), branchNames.end(), "time") != branchNames.end())
        tree->SetBranchAddress("time", &dataTime);
    else
        Fatal("collectData", "Requested branch 'time' does not exist in the file.");

    // Get the number of entries
    int entries = tree->GetEntries();

    // Collect the data
    for (int i = 0; i < entries; i++) {
        // Get the corresponding entry
        tree->GetEntry(i);

        // Collect the data
        virtualdetectorIds.push_back(dataVirtualdetectorId);
        pdgIds.push_back(dataPdgId);
        energies.push_back(dataE);
        times.push_back(dataTime);
    };

    // Check if data has been collected
    if (virtualdetectorIds.size() == 0)
        Fatal("collectData", "No data was collected from this file");

    // Clear
    file->Close();
    delete file;
    std::cout << "Finished processing file " << fileName << std::endl;
    return;
};

void collectDetectorData(const std::string fileName, const std::string treeName, std::vector<double> &energies, std::vector<double> &times) {
    /*
        Description
            Collects all the required data from detector TTrees

        Arguments
            fileName - as documented in function "plot"
            treeName - as documented in function "plot"
            energies - vector of energies associated with the relevant dataset, in MeV
            times - vector of times associated with the relevant dataset, in ns

        Variables
            file - ROOT TFile interface
            branches - ROOT TFile interface to TTree branches
            branchNames - vector of branch names
            dataE - energy from file
            dataTime - time from file
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
    double dataE, dataTime;
    if (std::find(branchNames.begin(), branchNames.end(), "E") != branchNames.end()) // RETURNTOME - change "totE" to "E"
        tree->SetBranchAddress("E", &dataE);
    else
        Fatal("collectData", "Requested branch 'E' does not exist in the file.");
    if (std::find(branchNames.begin(), branchNames.end(), "time") != branchNames.end())
        tree->SetBranchAddress("time", &dataTime);
    else
        Fatal("collectData", "Requested branch 'time' does not exist in the file.");

    // Get the number of entries
    int entries = tree->GetEntries();

    // Collect the data
    for (int i = 0; i < entries; i++) {
        // Get the corresponding entry
        tree->GetEntry(i);

        // Collect the data
        energies.push_back(dataE);
        times.push_back(dataTime);
    };

    // Check if data has been collected
    if (energies.size() == 0)
        Fatal("collectData", "There has been no data collected, so no plots will be generated");

    // Clear
    file->Close();
    delete file;
    std::cout << "Finished processing file " << fileName << std::endl;
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

void makePlot(std::vector<double> electronEnergies, std::vector<double> electronTimes, std::vector<double> muonEnergies, std::vector<double> muonTimes, const bool shiftEMin, double eRed, double binWidthFull, double binWidthRed, double binWidth347, double binWidth844, double binWidth1809, const double signalAcceptance, const double scaleFactor, const bool highResolution, const bool timeCuts, const bool flashCut347, const bool convertMeVTokeV, const std::string title, const std::string particleName) {
    /*
        Description
            Generates both both the stacked and unstacked plots for the full energy range, reduced energy range, and signal range

        Arguments
            electronEnergies - documented in function "plot"
            electronTimes - documented in function "plot"
            muonEnergies - documented in function "plot"
            muonTimes - documented in function "plot"
            shiftEMin - if true, sets the plot minimum to the first bin. Useful with the detector plots as many events have zero energy deposit, making the plots not display sufficient information
            eRed - documented in function "plot"
            binWidthFull - documented in function "plot"
            binWidthRed - documented in function "plot"
            binWidth347 - documented in function "plot"
            binWidth844 - documented in function "plot"
            binWidth1809 - documented in function "plot"
            signalAcceptance - documented in function "plot"
            scaleFactor - documented in function "plot"
            highResolution - documented in function "makePlots"
            timeCuts - documented in function "makePlots"
            flashCut347 - documented in function "makePlots"
            convertMeVTokeV - documented in function "makePlots"
            title - documented in function "makePlots"
            particleName - documented in function "makePlots"

        Variables
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
            vars - temporary buffer to store data if it requires unit conversion from MeV to keV
            conversion - conversion factor from MeV to keV
            eMaxElectron - maximum energy of data from EleBeamCat
            eMaxMuon - maximum energy of data from MuBeamCat
            nBinsFull - number of bins in full energy plot
            nBinsRed - numeber of bins in reduced energy plot
            nBins347 - number of bins in 347keV plot
            nBins844 - number of bins in 844keV plot
            nBins1809 - number of bins in 1809keV plot
            binWidthFullStr - bin width of full energy plot as a string
            binWidthRedStr - bin width of reduced energy plot as a string
            binWidth347Str - bin width of 347keV signal plot as a string
            binWidth844Str - bin width of 844keV signal plot as a string
            binWidth1809Str - bin width of 1809keV signal plot as a string
            eMinFull - minimum energy plotted in full energy plot
            eMinRed - minimum energy plotted in reduced energy plot
            order - order of plot generation histograms
            eRanges - energy ranges of histograms
            nBins - bin numbers in histograms
            binWidthStrs - bin widths as strings for histograms
            nOrder - number of plot categories
            i - iterator variable
            unit - energy unit as a string
            plotFileNamePrefix - prefix for plot file names
            plotFileNameSuffix - suffix for plot file names
            plotFileNameSuffixSignal - suffix for signal plot file names
            plotFileNameExtension - encoding used for plot files
            uFilenames - vector of file names of unstacked plots
            sFileNames - vector of file names of stacked plots
            plotTitle - vector of title string for plots
            uTitles - vector of titles of unstacked plots
            sTitles - vector of titles of stacked plots
            px - number of x pixels for canvases
            py - number of y pixels for canvases
            uCanvases - vector of canvases for unstacked plots
            sCanvases - vector of canvases for stacked plots
            Canvases - vector of all canvases
            lowResAxisTextSize - low resolution axis text size
            c - canvas iterator
            uHists - vector of unstacked histograms
            sHists - vector of stacked histograms
            nElectronPoints - number of data points from EleBeamCat
            e - energy buffer vector
            t - time buffer vector
            appliedScaling - monitors whether a scale factor was applied or not
            nMuonPoints - number of data points from MuBeamCat
            count - buffer vector for histogram bin checks
            lx1 - defines an x coordinate for stat boxes and legend
            ly1 - defines a y coordinate for stat boxes and legend
            lx2 - defines the other x coordinate for stat boxes and legend
            ly2 - defines the other y coordinate for stat boxes and legend
            legend - defines the TLegend
            o - `order` iterator
            hstacks - vector of all the THStacks
    */
    // Perform checks
    if (electronEnergies.empty() && muonEnergies.empty()) {
        std::cout << "Passed energy vectors are both empty, not generating histograms." << std::endl;
        return;
    };
    if ((flashCut347) && (!timeCuts))
        Fatal("plot", "If the flash cut is applied, the time cuts need to be applied too (if flashCut347=true, timeCuts=true is required)");
    if ((scaleFactor - 1) > std::numeric_limits<double>::epsilon() && title.substr(0, 2) != "VD")
        Fatal("makePlot", "If scaling is applied, detector data cannot be used");

    std::cout << std::string("Generating plots for ") + (particleName != "all" ? particleName : "all particles") + " at " + title + " with "  + (highResolution ? "high resolution, " : "low resolution, ") + (timeCuts ? "time cuts, " : "no time cuts, ") + (convertMeVTokeV ? "in keV" : "in MeV") << std::endl;

    // Define the signal parameters, energy in MeV and time in ns
    double e347  = 0.347,   eRange347  = e347  * signalAcceptance,  eMin347  = e347  - eRange347,   eMax347  = e347  + eRange347;
    double e844  = 0.844,   eRange844  = e844  * signalAcceptance,  eMin844  = e844  - eRange844,   eMax844  = e844  + eRange844;
    double e1809 = 1.809,   eRange1809 = e1809 * signalAcceptance,  eMin1809 = e1809 - eRange1809,  eMax1809 = e1809 + eRange1809;

    // Double the energy range to include both sides of the acceptance
    eRange347  *=2;
    eRange844  *=2;
    eRange1809 *=2;

    // Define the time ranges for the signals, in ns
    const double tMicroSpill = 1695;
    const double tMin347  = flashCut347 ? 200 : 300,    tMax347  = 700,     tMod347  = tMicroSpill;
    const double tMin844  = 492000,                     tMax844  = 1330000, tMod844  = tMax844;
    const double tMin1809 = 500,                        tMax1809 = 1600,    tMod1809 = tMicroSpill;

    // Convert everything to keV if needed
    double* vars[] = {&e347,  &eMin347,  &eMax347,  &eRange347,  &binWidth347,
                      &e844,  &eMin844,  &eMax844,  &eRange844,  &binWidth844,
                      &e1809, &eMin1809, &eMax1809, &eRange1809, &binWidth1809,
                      &eRed, &binWidthRed, &binWidthFull};
    if (convertMeVTokeV) {
        double conversion = 1000;
        std::transform(std::begin(vars),         std::end(vars),         std::begin(vars),         [](double* x)-> double*{ *x *= 1000; return x;});
        std::transform(electronEnergies.begin(), electronEnergies.end(), electronEnergies.begin(), [conversion](double x) { return x * conversion; });
        std::transform(muonEnergies.begin(),     muonEnergies.end(),     muonEnergies.begin(),     [conversion](double x) { return x * conversion; });
    };

    // Setup histograms bin counts
    double eMaxElectron = (electronEnergies.empty() ?   eRed : *std::max_element(electronEnergies.begin(), electronEnergies.end()));
    double eMaxMuon     = (muonEnergies.empty() ?       eRed : *std::max_element(muonEnergies.begin(),     muonEnergies.end()));
    double eMax         =  std::max(eMaxElectron, eMaxMuon);
    int nBinsFull = eMax/binWidthFull,
        nBinsRed  = eRed/binWidthRed,
        nBins347  = eRange347/binWidth347,
        nBins844  = eRange844/binWidth844,
        nBins1809 = eRange1809/binWidth1809;
    std::string binWidthFullStr = convertBinWidthToStr(binWidthFull),
                binWidthRedStr  = convertBinWidthToStr(binWidthRed),
                binWidth347Str  = convertBinWidthToStr(binWidth347),
                binWidth844Str  = convertBinWidthToStr(binWidth844),
                binWidth1809Str = convertBinWidthToStr(binWidth1809);

    // Define the spectra minima, offset useful for detector as zero entry strongly dominates
    const double eMinFull = (shiftEMin ? binWidthFull : 0);
    const double eMinRed  = (shiftEMin ? binWidthRed  : 0);

    // Set up the variables for setting up the histograms
    std::vector<std::string> order              = {"Full",              "Red",              "347",              "844",              "1809"};                // Order that TCanvases, TH1Ds, and THStacks are generated
    std::vector<std::vector<double>> eRanges    = {{eMinFull, eMax},    {eMinRed, eRed},    {eMin347, eMax347}, {eMin844, eMax844}, {eMin1809, eMax1809}};  // Energy ranges for TH1Ds
    std::vector<int> nBins                      = {nBinsFull,           nBinsRed,           nBins347,           nBins844,           nBins1809};             // Bin count for TH1Ds
    std::vector<std::string> binWidthStrs       = {binWidthFullStr,     binWidthRedStr,     binWidth347Str,     binWidth844Str,     binWidth1809Str};       // Bin widths as strings for plot titles
    const int nOrder = order.size(), sigRange = 1;

    // Validate the bin counts
    for (int i = 0; i < nOrder; i++) {
        if (nBins[i] < 2) {
            std::cout << "For plot '" << order[i] << "', energy range: ()" << eRanges[i][0] << ", " << eRanges[i][1] << "), binWidth" << binWidthStrs[i] << ", bin count" << nBins[i] << std::endl;
            Fatal("makePlot", "For the details above, it is forbidden to have a bin count of 1");
        };
    };

    // Set up the file names
    std::string unit = (convertMeVTokeV ? "keV" : "MeV");
    std::string plotFileNamePrefix = title + "." + particleName + ".";
    if ((title.substr(0, 2) == "VD") && ((scaleFactor - 1) > std::numeric_limits<double>::epsilon()))
        plotFileNamePrefix += "scaled.";
    std::string plotFileNameSuffix          = std::string("") + unit + "." + (highResolution ? "high" : "low");
    std::string plotFileNameSuffixSignal    = std::string("") + (timeCuts ? ".time-cut" : ".no-time-cut" ) + (flashCut347 ? ".347-flash" : "");
    std::string plotFileNameExtension       = ".png";
    // Store the file names
    std::vector<std::string> uFileNames, sFileNames;
    for (int i = 0; i < nOrder; i++) {
        uFileNames.emplace_back(plotFileNamePrefix + order[i] + ".no-split." + plotFileNameSuffix + (i > sigRange ? plotFileNameSuffixSignal : "" ) + plotFileNameExtension);
        sFileNames.emplace_back(plotFileNamePrefix + order[i] + ".split."    + plotFileNameSuffix + (i > sigRange ? plotFileNameSuffixSignal : "" ) + plotFileNameExtension);
    };

    // Set up the plot title
    std::string plotTitle = "";
    if (std::abs(scaleFactor - 1) > std::numeric_limits<double>::epsilon())
        plotTitle = "Scaled " + particleName;
    else
        plotTitle = (char)std::toupper(particleName[0]) + particleName.substr(1);
    plotTitle += (particleName == "all" ? " particles at " : "s at ") + title + (title.substr(0, 2) != "VD" ? " detector" : "");
    // store the plot titles
    std::vector<std::string> uTitles, sTitles;
    std::string titleStr = "";
    for (int i = 0; i < nOrder; i++) {
        titleStr = plotTitle + "; Kinetic Energy [" + unit + "];Count / " + binWidthStrs[i] + " " + unit;
        uTitles.emplace_back(titleStr);
        sTitles.emplace_back(titleStr);
    };

    // Set up the TCanvases
    int px = (highResolution ? 1500 : 750), py = (highResolution ? 1000 : 500);
    std::vector<TCanvas*> uCanvases, sCanvases;
    for (std::string o : order) {
        uCanvases.emplace_back(new TCanvas(("c" +o).c_str(), ("c" +o).c_str(), px, py)); // Unstacked TCanvas
        sCanvases.emplace_back(new TCanvas(("cs"+o).c_str(), ("cs"+o).c_str(), px, py)); // Stacked TCanvas
    };
    // Store the canvases
    std::vector<TCanvas*> canvases(uCanvases.begin(), uCanvases.end());
    canvases.insert(canvases.end(), sCanvases.begin(), sCanvases.end());

    // Apply TCanvas formatting
    double lowResAxisTextSize = 0.06;
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

    // Set up the TH1Ds
    std::vector<TH1D*> uHists;
    std::vector<std::vector<TH1D*>> sHists;
    for (int i = 0; i < nOrder; i++) {
        uHists.emplace_back(
            new TH1D(
                ("Combined" + order[i]).c_str(),
                ("Combined" + order[i]).c_str(),
                nBins[i], eRanges[i][0], eRanges[i][1])
        );
        sHists.emplace_back(std::vector<TH1D*>{
            new TH1D(
                ("Ele" + order[i]).c_str(),
                ("Ele" + order[i]).c_str(),
                nBins[i], eRanges[i][0], eRanges[i][1]),
            new TH1D(
                ("Mu" + order[i]).c_str(),
                ("Mu" + order[i]).c_str(),
                nBins[i], eRanges[i][0], eRanges[i][1])
        });
    };
    for (std::vector<TH1D*> h : sHists) {
        for (TH1D* sh : h)
            sh->SetStats(0);
    };

    // Populate the histograms with electron data
    int nElectronPoints = electronEnergies.size();
    double e = 0.0, t = 0.0;
    for (int i = 0; i < nElectronPoints; i++) {
        t = electronTimes[i];
        e = electronEnergies[i];
        if (eMinFull < e && e < eMax) {
            uHists[0]->Fill(e);
            sHists[0][0]->Fill(e);
        };
        if (eMinRed < e  && e < eRed) {
            uHists[1]->Fill(e);
            sHists[1][0]->Fill(e);
        };
        if ((!timeCuts && eMin347 < e  && e < eMax347)  || (timeCuts && tMin347 < fmod(t, tMod347)   && fmod(t, tMod347) < tMax347   && eMin347 < e  && e < eMax347))  {
            uHists[2]->Fill(e);
            sHists[2][0]->Fill(e);
        };
        if ((!timeCuts && eMin844 < e  && e < eMax844)  || (timeCuts && tMin844 < fmod(t, tMod844)   && fmod(t, tMod844) < tMax844   && eMin844 < e  && e < eMax844))  {
            uHists[3]->Fill(e);
            sHists[3][0]->Fill(e);
        };
        if ((!timeCuts && eMin1809 < e && e < eMax1809) || (timeCuts && tMin1809 < fmod(t, tMod1809) && fmod(t, tMod1809) < tMax1809 && eMin1809 < e && e < eMax1809)) {
            uHists[4]->Fill(e);
            sHists[4][0]->Fill(e);
        };
    };

    // Scale the electron data if appropriate
    bool appliedScaling = false;
    if ((scaleFactor - 1) > std::numeric_limits<double>::epsilon()) {
        appliedScaling = true;
        for (TH1D* h : uHists)
            h->Scale(scaleFactor);
        for (std::vector<TH1D*> h : sHists)
            h[0]->Scale(scaleFactor);
    };

    // Populate the histograms with muon data
    int nMuonPoints = muonEnergies.size();
    for (int i = 0; i < nMuonPoints; i++) {
        t = muonTimes[i];
        e = muonEnergies[i];
        if (eMinFull < e && e < eMax) {
            uHists[0]->Fill(e);
            sHists[0][1]->Fill(e);
        };
        if (eMinRed < e && e < eRed) {
            uHists[1]->Fill(e);
            sHists[1][1]->Fill(e);
        };
        if ((!timeCuts && eMin347 < e  && e < eMax347)  || (timeCuts && tMin347 < fmod(t, tMod347)   && fmod(t, tMod347) < tMax347   && eMin347 < e  && e < eMax347))  {
            uHists[2]->Fill(e);
            sHists[2][1]->Fill(e);
        };
        if ((!timeCuts && eMin844 < e  && e < eMax844)  || (timeCuts && tMin844 < fmod(t, tMod844)   && fmod(t, tMod844) < tMax844   && eMin844 < e  && e < eMax844))  {
            uHists[3]->Fill(e);
            sHists[3][1]->Fill(e);
        };
        if ((!timeCuts && eMin1809 < e && e < eMax1809) || (timeCuts && tMin1809 < fmod(t, tMod1809) && fmod(t, tMod1809) < tMax1809 && eMin1809 < e && e < eMax1809)) {
            uHists[4]->Fill(e);
            sHists[4][1]->Fill(e);
        };
    };

    // Check the hists for overflow and underflow bins
    int count = 0;
    for (int i = 0; i < nOrder; i++) {
        count = uHists[i]->GetBinContent(0);
        if (count) {
            std::cout << "In plot " << uFileNames[i] << ", number of underflow bins is " << count << std::endl;
            std::cout << "This histogram's lower range is: " << uHists[i]->GetXaxis()->GetXmin() << std::endl;
            Fatal("makePlot", "Generated plot has values in the underflow bins");
        };
        count = uHists[i]->GetBinContent(uHists[i]->GetNbinsX() + 1);
        if (count) {
            std::cout << "In plot " << uFileNames[i] << ", number of overflow bins is " << count << std::endl;
            std::cout << "This histogram's upper range is: " << uHists[i]->GetXaxis()->GetXmax() << std::endl;
            Fatal("makePlot", "Generated plot has values in the overflow bins");
        };
        count = uHists[i]->GetEntries() - uHists[i]->Integral(1, uHists[i]->GetNbinsX());
        if (!appliedScaling && count) { // scaling affects integral but not entries
            std::cout << "In plot " << uFileNames[i] << ", difference between entries and integral is " << count << std::endl;
            Fatal("makePlot", "Generated plot has values in the underflow/overflow bins");
        };
        count = sHists[i][0]->GetBinContent(0) + sHists[i][1]->GetBinContent(0);
        if (count) {
            std::cout << "In plot " << sFileNames[i] << ", number of underflow bins is " << count << std::endl;
            std::cout << "This histogram's lower range is: " <<sHists[i][0]->GetXaxis()->GetXmin() << std::endl;
            Fatal("makePlot", "Generated plot has values in the underflow bins");
        };
        count = sHists[i][0]->GetBinContent(sHists[i][0]->GetNbinsX() + 1) + sHists[i][1]->GetBinContent(sHists[i][1]->GetNbinsX() + 1);
        if (count) {
            std::cout << "In plot " << sFileNames[i] << ", number of overflow bins is " << count << std::endl;
            std::cout << "This histogram's upper range is: " << sHists[i][0]->GetXaxis()->GetXmax() << std::endl;
            Fatal("makePlot", "Generated plot has values in the overflow bins");
        };
        count = sHists[i][0]->GetEntries() - sHists[i][0]->Integral(1, sHists[i][0]->GetNbinsX()) + sHists[i][1]->GetEntries() - sHists[i][1]->Integral(1, sHists[i][1]->GetNbinsX());
        if (!appliedScaling && count) { // scaling affects integral but not entries
            std::cout << "In plot " << sFileNames[i] << ", difference between entries (" << sHists[i][0]->GetEntries() + sHists[i][1]->GetEntries() << ") and integral ( " << sHists[i][0]->Integral(1, sHists[i][0]->GetNbinsX()) + sHists[i][1]->Integral(1, sHists[i][1]->GetNbinsX()) << ") is " << count << std::endl;
            Fatal("makePlot", "Generated plot has values in the underflow/overflow bins");
        };
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

    // Draw and save the unstacked plots
    for (int i = 0; i < nOrder; i++) {
        uCanvases[i]->cd();
        count = uHists[i]->GetEntries();
        if (count == 0) {
            std::cout << uFileNames[i] << " has no entries, not saving" << std::endl;
            continue;
        };
        if (appliedScaling) // When scaling, the number of entries is not scaled
            uHists[i]->SetEntries(uHists[i]->Integral(1, uHists[i]->GetNbinsX()));
        uHists[i]->SetTitle(uTitles[i].c_str());
        uHists[i]->SetMinimum(0); // Sets the y minimum
        uHists[i]->Draw("HIST");
        if (!highResolution) {
            uHists[i]->GetXaxis()->SetLabelSize(lowResAxisTextSize);
            uHists[i]->GetXaxis()->SetTitleSize(lowResAxisTextSize);
            uHists[i]->GetYaxis()->SetLabelSize(lowResAxisTextSize);
            uHists[i]->GetYaxis()->SetTitleSize(lowResAxisTextSize);
        };
        uCanvases[i]->SaveAs(uFileNames[i].c_str());
    };
    // Clean up from the unstacked plots
    for (int i = nOrder - 1; i > -1; i--) {
        uHists[i]->Delete();
        uCanvases[i]->Close();
    };

    // Set up the legend
    TLegend *legend = new TLegend(lx1, ly1, lx2, ly2, "Data");
    legend->SetHeader("Dataset", "C");
    legend->AddEntry("Background",  "Background",   "l");
    legend->AddEntry("Signal",      "Signal",       "l");
    for (TCanvas* c : sCanvases)
        c->Update();

    // Set up the THStacks
    std::vector<THStack*> hStacks;
    for (std::string o : order)
        hStacks.emplace_back(new THStack(("hs" + o).c_str(), ("hs" + o).c_str()));

    // Save the stacked plots
    for (int i = 0; i < nOrder; i++) {
        sCanvases[i]->cd();
        count = sHists[i][0]->GetEntries() + sHists[i][1]->GetEntries();
        if (count == 0) {
            std::cout << "Plot " << sFileNames[i] << " has no entries, not generating" << std::endl;
            continue;
        };
        hStacks[i]->Add(sHists[i][0]);
        hStacks[i]->Add(sHists[i][1]);
        sHists[i][0]->SetFillColor(kRed);
        sHists[i][0]->SetLineColor(kRed);
        sHists[i][1]->SetFillColor(kBlue);
        sHists[i][1]->SetLineColor(kBlue);
        hStacks[i]->SetTitle(sTitles[i].c_str());
        hStacks[i]->SetMinimum(0); // Sets the y minimum
        hStacks[i]->Draw("HIST");
        legend->Draw();
        if (!highResolution) {
            hStacks[i]->GetXaxis()->SetLabelSize(lowResAxisTextSize);
            hStacks[i]->GetXaxis()->SetTitleSize(lowResAxisTextSize);
            hStacks[i]->GetYaxis()->SetLabelSize(lowResAxisTextSize);
            hStacks[i]->GetYaxis()->SetTitleSize(lowResAxisTextSize);
        };
        sCanvases[i]->SaveAs(sFileNames[i].c_str());
    };

    // Cleanup
    for (int i = nOrder - 1; i > -1; i--) {
        sHists[i][1]->Delete();
        sHists[i][0]->Delete();
        hStacks[i]->Delete();
        sCanvases[i]->Close();
    };

    std::cout << "Finished\n" << std::endl;
    return;
};

void makePlots(std::vector<double> &electronEnergies, std::vector<double> &electronTimes, std::vector<double> &muonEnergies, std::vector<double> &muonTimes, const bool shiftEMin, double eRed, double binWidthFull, double binWidthRed, double binWidth347, double binWidth844, double binWidth1809, const double signalAcceptance, const double scaleFactor, const std::string title, const std::string particleName) {
    /*
        Description
            Generates all the required plots

        Arguments - all parameters passed as references to conserve memory
            electronEnergies - as documented in function `plot`
            electronTimes - as documented in function `plot`
            muonEnergies - as documented in function `plot`
            muonTimes - as documented in function `plot`
            eRed - as documented in function `plot`
            binWidthFull - as documented in function `plot`
            binWidthRed - as documented in function `plot`
            title - prefix for the title, this will either be "VD" for virtual detectors or either "HPGe" or "LaBr" for the detectors. This is also used for the plot file name
            particleName - name of the particle being plotted, must be one of plotParticleNames or "all". If a different particle is required to be plotted, add it
            scaleFactor - as documented in function `plot`

        Variables
            plotTitles - available plot titles
            plotTitle - iterator for `plotTitles`
            boolValues - vector of valid boolean values
            scaleFactors - values of scale factors to apply with function `makePlot`
            convertMevTokeV - variable that controls whether the plot should be generated in MeV or keV
            scale - value of scale factor to apply with function `makePlot`
            highRes - variable that controls whether the plot should be generated in high or low resolution
    */
    // Perform pre plot checks on the plot title and file name
    std::vector<std::string> plotTitles = {"VD", "HPGe", "LaBr"};
    if ((std::find(plotTitles.begin(), plotTitles.end(), title) == plotTitles.end()) && (title.substr(0, 2) != "VD")) {
        std::cout << "The plot requested with title " << title << " is not one of ";
        for (std::string plotTitle : plotTitles)
            std::cout << plotTitle << ", ";
        std::cout << "not generating this plot." << std::endl;
        Fatal("makePlots", "Incorrect plot title format!");
    };

    // Construct loop variables for plotting
    std::vector<bool> boolValues = {true, false};
    std::vector<double> scaleFactors = {1.0};
    if (std::abs(scaleFactor - 1) > std::numeric_limits<double>::epsilon())
        scaleFactors.emplace_back(scaleFactor);

    // Call the plotting function
    for (bool convertMeVTokeV : boolValues) {
        for (double scale : scaleFactors) {
            for (bool highRes : boolValues) {
                makePlot(electronEnergies, electronTimes, muonEnergies, muonTimes, shiftEMin, eRed, binWidthFull, binWidthRed, binWidth347, binWidth844, binWidth1809, signalAcceptance, scale, highRes,  false, false, convertMeVTokeV, title, particleName);
                makePlot(electronEnergies, electronTimes, muonEnergies, muonTimes, shiftEMin, eRed, binWidthFull, binWidthRed, binWidth347, binWidth844, binWidth1809, signalAcceptance, scale, highRes,  true,  false, convertMeVTokeV, title, particleName);
                makePlot(electronEnergies, electronTimes, muonEnergies, muonTimes, shiftEMin, eRed, binWidthFull, binWidthRed, binWidth347, binWidth844, binWidth1809, signalAcceptance, scale, highRes,  true,  true , convertMeVTokeV, title, particleName);
            };
        };
    };

    return;
};

void plotVirtualdetector(std::vector<double> &electronEnergies, std::vector<double> &electronTimes, const std::vector<ULong64_t> &electronVirtualdetectorIds, const std::vector<int> &electronPdgIds, std::vector<double> &muonEnergies, std::vector<double> &muonTimes, const std::vector<ULong64_t> &muonVirtualdetectorIds, const std::vector<int> &muonPdgIds, double eRed, double binWidthFull, double binWidthRed, double binWidth347, double binWidth844, double binWidth1809, const double signalAcceptance, const double scaleFactor) {
    /*
        Description
            Generates all the plots for virtual detectors

        Arguments
            electronEnergies - as documented in function `plot`
            electronTimes - as documented in function `plot`
            electronVirtualdetectorIds - as documented in function `plot`
            electronPdgIds - as documented in function `plot`
            muonEnergies - as documented in function `plot`
            muonTimes - as documented in function `plot`
            muonVirtualdetectorIds - as documented in function `plot`
            muonPdgIds - as documented in function `plot`
            eRed - as documented in function `plot`
            binWidthFull - as documented in function `plot`
            binWidthRed - as documented in function `plot`

        Variables
            plotVirtualDetectorIds - virtual detector IDs for which plots will be generated
            plotPdgIds - PDG IDs for which plots will be generated
            nPlotPdgIds - number of PDG IDs to plot for
            plotPdgId - PDG ID to plot
            plotTitle - plot title start
            plotParticleNames - names of particles to plot for
            plotParticleName - name of particle being plotted
            plotElectronEnergies - vector of electron energies to plot
            plotElectronTimes - vector of electron times to plot
            plotMuonEnergies - vector of muon energies to plot
            plotMuonTimes - vector of muon times to plot
            nElectronPoints - number of data points from EleBeamCat data
            nMuonPoints - number of data points from MuBeamCat data
            timeCuts - if true, applies a modulus operator to the time and applies the standard time cuts of 300ns < t_347 < 700ns, 492ms < t_844 < 1330ms, 500ns < t_1809 < 1600ns. If false, no time cuts are applied
            flashCut347 - if true, applies a 200ns < t_347 < 700ns time cut instead of 300ns < t_347 < 700ns. This requires timeCuts to be true
            eRed - maximum energy of reduced spectrum to plot
            virtualdetectorId - selects a single virtual detector to plot. If left to 0, will plot all the virtual detectors available
            eRed - sets the maximum energy to plot for the reduced spectrum in MeV
            plotKeV - if true generates plots in keV instead of MeV
    */

    // Make a unique list of virtual detector IDs
    // std::set<ULong64_t> plotVirtualdetectorIds(electronVirtualdetectorIds.begin(), electronVirtualdetectorIds.end());
    // for (ULong64_t plotVirtualdetectorId : muonVirtualdetectorIds)
    //     plotVirtualdetectorIds.emplace(plotVirtualdetectorId);
    // std::cout << "Found virtualdetectors with ID: ";
    // for (double plotVirtualdetectorId : plotVirtualdetectorIds)
    //     std::cout << plotVirtualdetectorId << ", ";
    // std::cout << "\n" << std::endl;
    std::set<ULong64_t> plotVirtualdetectorIds = {88, 89, 90, 101};

    // Make a unique list of PDG IDs
    // std::set<int> plotPdgIds(electronPdgIds.begin(), electronPdgIds.end());
    // for (int plotPdgId : muonPdgIds)
    //     plotPdgIds.emplace(plotPdgId);
    // std::cout << "Found particles with PDG ID: ";
    // for (int plotPdgId : plotPdgIds)
    //     std::cout << plotPdgId << ", ";
    // std::cout << "\n" << std::endl;
    std::vector<int> plotPdgIds = {0, -11, 11, 22, 2112};
    int nPlotPdgIds = plotPdgIds.size(), plotPdgId = 0;
    std::vector<std::string> plotParticleNames = {"all", "positron", "electron", "photon", "neutron"};
    std::string plotParticleName = "";

    // Set up plotting parameters
    std::string plotTitle = "";
    std::vector<double> plotElectronEnergies, plotElectronTimes, plotMuonEnergies, plotMuonTimes;

    // Make plot by Virtual Detector ID
    int nElectronPoints = electronEnergies.size(), nMuonPoints = muonEnergies.size();
    for (ULong64_t plotVirtualdetectorId : plotVirtualdetectorIds) {
        plotTitle = "VD" + std::to_string(plotVirtualdetectorId);
        // Make plot by PDG ID
        for (int i = 0; i < nPlotPdgIds; i++) {
            plotPdgId = plotPdgIds[i];
            plotParticleName = plotParticleNames[i];
            std::cout << "Collecting data for virtualdetector ID " << plotVirtualdetectorId << " and particle " << plotParticleName << std::endl;

            // Select the appropriate dat for the plots
            for (int j = 0; j < nElectronPoints; j++) {
                if (electronVirtualdetectorIds[j] == plotVirtualdetectorId && (plotPdgId == 0 || electronPdgIds[j] == plotPdgId)) {
                    plotElectronEnergies.emplace_back(electronEnergies[j]);
                    plotElectronTimes.emplace_back(electronTimes[j]);
                };
            };
            for (int j = 0; j < nMuonPoints;     j++) {
                if (muonVirtualdetectorIds[j] == plotVirtualdetectorId     && (plotPdgId == 0 || muonPdgIds[j] == plotPdgId)) {
                    plotMuonEnergies.emplace_back(muonEnergies[j]);
                    plotMuonTimes.emplace_back(muonTimes[j]);
                };
            };
            makePlots(plotElectronEnergies, plotElectronTimes, plotMuonEnergies, plotMuonTimes, false, eRed, binWidthFull, binWidthRed, binWidth347, binWidth844, binWidth1809, signalAcceptance, scaleFactor, plotTitle, plotParticleName);
            plotElectronEnergies.clear();
            plotElectronTimes.clear();
            plotMuonEnergies.clear();
            plotMuonTimes.clear();
        };
    };
    return;
};

void plotDetector(std::vector<double> &electronEnergies, std::vector<double> &electronTimes, std::vector<double> &muonEnergies, std::vector<double> &muonTimes, std::string detectorName, double eRed, double binWidthFull, double binWidthRed, double binWidth347, double binWidth844, double binWidth1809, const double signalAcceptance) {
    /*
        Description
            Generates plots for the detectors

        Arguments
            ekectronEnergies - as documented in function `plot`
            electronTimes - as documented in function `plot`
            muonEnergies - as documented in function `plot`
            muonTimes - as documented in function `plot`
            detectorName - as documented in function `plot`
            eRed - as documented in function `plot`
            binWidthFull - as documented in function `plot`
            binWidthRed - as documented in function `plot`
            binWidth347 - as documented in function `plot`
            binWidth844 - as documented in function `plot`
            binWidth1809 - as documented in function `plot`
            signalAcceptance - as documented in function `plot`
    */

    makePlots(electronEnergies, electronTimes, muonEnergies, muonTimes, true, eRed, binWidthFull, binWidthRed, binWidth347, binWidth844, binWidth1809, signalAcceptance, 1, detectorName, "all");
    return;
};

void plotAllSpectra(const std::vector<std::string> electronFileNames, const std::vector<std::string> muonFileNames, const std::string treeName, const std::string virtualdetectorOrDetectorName, const double scaleFactor = 1, const double signalAcceptance = 0.1, double eRed = 2.0, double binWidthFull = 0.5, double binWidthRed = 0.01, double binWidth347 = 0.005, double binWidth844 = 0.005, double binWidth1809  = 0.01) {
    /*
        Description
            Plots the spectra from STM simulation campaigns for electrons, positrons, photons, and neutrons. For other particles, the associated PDG IDs and particle names need to be added to the vectors plotPdgIds and plotParticleNames
            Uses the results of MakeTree.fcl
            usage example - see plot.sh

        Notes
            Plot file names
              VD<virtual detector ID>.<particle name>.<scaled/>.<Full/Red/signal>.<split/no-split>.<keV/MeV>.<high/low>.<time-cut/no-time-cut>.<347-flash/>.png
              <detector>.<particle.name>.<full/red/signal>.<split/no-split>.<keV/MeV>.<high/low>.<time-cut/no-time-cut>.<347-flash/>.png
            Plot titles
              <scaled/> <particle name>s  at VD<virtual detctor ID>
              <particle name>s at <detector> detector
            Plot file name and title parameters
                <virtual detector ID> - ID of the virtual detector
                <particle name> - name of the particle as defined in plotPdgIds and plotParticleNames, or "all" if all particles are used
                <scaled/> - if "scaled", scales the EleBeamCat results up by scaleFactor to the same POT statistics as the MuBeamCat results, if "", no scaling is applied. The scaling value is defined in plotScaling, and should only be used with the Stage 1 results
                <Full/Red/signal> - if "Full", the full energy range is plotted, and if "Red", the reduced energy range is plotted up to eRed, if "signal" the signal range is plotted ("signal" will be one of {347, 844, 1809})
                <split/no-split> - if "split", the spectrum is plotted splitting the results from EleBeamCat and MuBeamCat, and if "no-split", a single histogram is generated
                <keV/MeV> - if "keV", the spectrum is plotted in keV, if "MeV", the spectrum is plotted in MeV
                <high/low> - if "high", the spectrum is plotted in high resolution, if "low", the spectrum is plotted in low resolution
                <time-cut/no-time-cut> - if "time-cut", the time cut of 300ns < t_347 < 700ns, 492ms < t_844 < 1330ms, 500ns < t_1809 < 1600ns is applied, and if "no-time-cut", no time cut is applied
                <347-flash/> - if "347-flash", the time cut of 200ns < t_347 < 700ns
                <detector> - name of the detector

            Required hard-coded changes for the user
            For virtual detectors, if plots are required for virtual detectors other than {101, 88, 89, 90}, this can be changed with variable plotVirtualdetectorIds
            For all plots, if plots are required for particles other than {positrons, electrons, muons, photons}, this can be changed with BOTH plotPdgIds and plotParticleNames

        Arguments
            electronFileNames - vector of EleBeamCat derived ROOT file names as a relative path to cwd
            muonFileNames - vector of MuBeamCat derived ROOT file names as a relative path to cwd
            treeName - name of tree containing all the data for electronFileNames and muonFileNames
            virtualdetectorOrDetector - string as either "virtualdetector" or the name of a STM detector, either "HPGe" or "LaBr". Selects the data to plot and plots it

        Optional arguments
            scaleFactor - scales the EleBeamCat results to make the EleBeamCat and MuBeamCat POT statistic equal by this amount
            signalAcceptance - amount of the energy window to include in the signal plots, default = 0.1 = 10%
            eRed - maximum energy of the reduced spectrum, in MeV
            binWidthFull - bin width of the full spectra, in MeV
            binWidthRed - bin width of the reduced spectra, in MeV
            binWidth347 - width of the bins in the 347keV spectrum in MeV
            binWidth844 - width of the bins in the 844keV spectrum in MeV
            binWidth1809 - width of the bins in the 1809keV spectrum in MeV

        Variables
            electronsVirtualdetectorId - for electrons, vector of virtualdetector IDs
            electronPdgId - for electrons, set of PDG IDs
            electronEnergies - for electrons, vector of energies in MeV
            electronTimes - for electrons, vector of times in ns
            muonVirtualdetectorId - for muons, vector of virtualdetector IDs
            muonPdgId - for muons, set of PDG IDs
            muonEnergies - for muons, vector of energies in MeV
            muonTimes - for muons, vector of times in ns
            fileName - iterator for the file name vectors
    */

    // Update global parameters
    SetErrorHandler(customErrorHandler);
    gROOT->SetBatch(kTRUE);

    // Perform pre data handling checks
    std::vector<std::string> virtualdetectorOrDetectorNames{"virtualdetector", "HPGe", "LaBr"};
    if (std::find(virtualdetectorOrDetectorNames.begin(), virtualdetectorOrDetectorNames.end(), virtualdetectorOrDetectorName) == virtualdetectorOrDetectorNames.end()) {
        std::cout << "Provided virtualdetectorOrDetectorName: " << virtualdetectorOrDetectorName << std::endl;
        Fatal("plot", "virtualdetectorOrDetectorName has to be either 'virtualdetector' if plotting for virutaldetector data, or one of 'HPGe' or 'LaBr' for the STM detectors.");
    };
    if (virtualdetectorOrDetectorName != "virtualdetector" && (scaleFactor - 1) > std::numeric_limits<double>::epsilon())
        Fatal("plot", "If plotting for the detectors, expected to not use scaling");

    // Initialize the data collection variables
    std::vector<double> electronEnergies, electronTimes, muonEnergies, muonTimes;
    std::vector<ULong64_t> electronVirtualdetectorIds, muonVirtualdetectorIds;
    std::vector<int> electronPdgIds, muonPdgIds;

    // Collect the data for each file called fileName (each file in the vector of electronFileNames and muonFileNames)
    std::cout << "Collecting data" << std::endl;
    if (virtualdetectorOrDetectorName == "virtualdetector") {
        for (const std::string fileName : electronFileNames)
            collectVirtualdetectorData(fileName, treeName, electronEnergies, electronTimes, electronVirtualdetectorIds, electronPdgIds);
        for (const std::string fileName : muonFileNames)
            collectVirtualdetectorData(fileName, treeName, muonEnergies,     muonTimes,     muonVirtualdetectorIds,     muonPdgIds);
        std::cout << "Data collected, generating plots\n" << std::endl;
        plotVirtualdetector(electronEnergies, electronTimes, electronVirtualdetectorIds, electronPdgIds, muonEnergies, muonTimes, muonVirtualdetectorIds, muonPdgIds, eRed, binWidthFull, binWidthRed, binWidth347, binWidth844, binWidth1809, signalAcceptance, scaleFactor);
    }
    else {
        for (std::string fileName : electronFileNames)
            collectDetectorData(fileName, treeName, electronEnergies, electronTimes);
        for (std::string fileName : muonFileNames)
            collectDetectorData(fileName, treeName, muonEnergies, muonTimes);
        std::cout << "Data collected, generating plots\n" << std::endl;
        plotDetector(electronEnergies, electronTimes, muonEnergies, muonTimes, virtualdetectorOrDetectorName, eRed, binWidthFull, binWidthRed, binWidth347, binWidth844, binWidth1809, signalAcceptance);
    };
   return;
};
