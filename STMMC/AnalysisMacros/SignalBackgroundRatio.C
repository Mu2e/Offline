// Calculates the signal to background ratio in the signal regions
// See SignalBackgroundRatio.sh for usage examples
// Original author - Pawel Plesniak

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
    if (!file || file->IsZombie())
        Fatal("collectData", "Failed to open the file.");

    // Get the tree
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


unsigned long count(std::vector<double> &energies, std::vector<double> &times, const double signalEnergy, const double signalAcceptance, std::vector<double> signalTimes) {
    /*
        Description
            Counts the number of entries in an enerrgy range with the relevant time cuts

        Variables
            energies - vector of energies to count within the energy windows
            times - vector of times associated with the given energy
            signalEnergy - energy of the signals [MeV]
            signalAcceptance - width of window to accept both signals and backgrounds in, as a multiple of the signal energy
            signalTimes - time cuts to apply to the signal as {tMin, tMax, tMod}. tMod corresponds to the modulo to apply to the signal time

        Variables
            eMin - minimum energy to accept [MeV]
            eMax - maximum energy to accept [MeV]
            tMin - minimum time to accept [ns]
            tMax - maximum time to accept [ns]
            tMod - modulus time to apply [ns]
            count - counter
            entries - number of entries in the energy and times vectors
            e - buffer variable for the energy
            t - buffer variable for the time
            tModded - time after the modulus has been applied
    */
    // Construct the limit variables
    const double eMin = signalEnergy * (1 - signalAcceptance), eMax = signalEnergy * (1 + signalAcceptance);
    const double tMin = signalTimes[0], tMax = signalTimes[1], tMod = signalTimes[2];
    unsigned long count = 0, entries = energies.size();

    // Construct the buffer variables
    double e = 0.0, t = 0.0, tModded = 0.0;

    // Count the relevant entries
    for (unsigned long i = 0; i < entries; i++) {
        e = energies[i];
        t = times[i];
        tModded = fmod(t, tMod);
        if (e > eMin && e < eMax && tModded > tMin && tModded < tMax)
            count++;
    };

    return count;
};

void SignalBackgroundRatio(const std::vector<std::string> electronFileNames, const std::vector<std::string> muonFileNames, const std::string treeName, const double signalAcceptance){
    /*
        Description
            Generates the signal to background ratio for each of the STM signal photons

        Arguments
            electronFileNames - vector of file names that generated the background to the signals
            muonFileNames - vector of file names that generated the signals
            treeName - name of the ttree containing the relevant data
            signalAcceptance - width of the window used to generate the signal to background ratio

        Variables
            electronEnergies - collected background energies
            electronTimes - collected background energy times
            muonEnergies - collected signal energies
            muonTimes - collected signal times
            fileName - iterator variable
            electronCount - count of background events in event window
            muonCount - count of signal events in event window
            signalEnergies - eneergies of STM signal photons
            signalTimes - vector of vectors of times used to accept the signals, as tMin, tMax, and tMod
            nSignals - number of signals to generate the signal to background ratio for
    */
    std::vector<double> electronEnergies, electronTimes, muonEnergies, muonTimes;
    for (std::string fileName : electronFileNames)
        collectDetectorData(fileName, treeName, electronEnergies, electronTimes);
    for (std::string fileName : muonFileNames)
        collectDetectorData(fileName, treeName, muonEnergies, muonTimes);

    // Initialize the data counter variables
    unsigned long electronCount = 0, muonCount = 0;
    std::vector<double> signalEnergies = {0.347, 0.844, 1.809};
    std::vector<std::vector<double>> signalTimes = {{300, 700, 1695}, {492000, 1330000, 1330000}, {500, 1600, 1695}}; // Arranged as tMin, tMax, tMod, s.t. tMod is the modulus applied to the time
    int nSignals = signalEnergies.size();
    std::cout << "Using a window acceptance of " << signalAcceptance << std::endl;
    for (int i = 0; i < nSignals; i++) {
        electronCount   = count(electronEnergies,   electronTimes,  signalEnergies[i],  signalAcceptance,   signalTimes[i]);
        muonCount       = count(muonEnergies,       muonTimes,      signalEnergies[i],  signalAcceptance,   signalTimes[i]);
        std::cout << "For signal at " << signalEnergies[i] << " MeV, the signal/background ratio is " << (1.0 * muonCount) / electronCount << std::endl;
    };
    return;
};
