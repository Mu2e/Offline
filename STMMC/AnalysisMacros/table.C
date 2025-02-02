// Generates summary tables of data generated from STM MC method, see doc db 00000 for more information
// See table.sh for usage examples
// Original author: Pawel Plesniak

void customErrorHandler(int level, Bool_t abort, const char* location, const char* message) {
    /*
        Description
        Define a custom error handler that won't print the stack trace but will print an error message and exit.
    */
    std::cerr << message << std::endl;
    if (level > kInfo)
        exit(1);
};

void collectVirtualdetectorData(const std::string &fileName, const std::string &treeName, const int &virtualdetectorId, std::vector<double> &energies, std::vector<int> &pdgIds) {
    /*
        Description
        Collects all the required data from virtual detector TTrees

        Arguments
        fileName - as documented in function "table"
        treeName - as documented in function "table"
        energies - vector of energies associated with the relevant dataset
        virtualdetectorId - as documented in function "table"
        pdgIds - vector of PDG IDs associated with the relevant dataset

        Variables
        file - ROOT TFile interface
        branches - ROOT TFile interface to TTree branches
        branchNames - vector of branch names
        dataVirtualDetectorId - virtual detector ID from file
        dataPdgId - PDG ID from file
        dataKE - kinetic energy from file
        entries - number of entries in the TTree
    */
    std::cout << "Processing file " << fileName << std::endl;

    // Get the branch
    TFile *file = new TFile(fileName.c_str());
    if (!file || file->IsZombie()) {
        Fatal("collectData", "Failed to open the file.");
    }
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
    double dataKE;
    if (std::find(branchNames.begin(), branchNames.end(), "virtualdetectorId") != branchNames.end())
        tree->SetBranchAddress("virtualdetectorId", &dataVirtualdetectorId);
    else
        Fatal("collectData", "Requested branch 'virtualdetectorId' does not exist in the file.");
    if (std::find(branchNames.begin(), branchNames.end(), "pdgId") != branchNames.end())
        tree->SetBranchAddress("pdgId", &dataPdgId);
    else
        Fatal("collectData", "Requested branch 'pdgId' does not exist in the file.");
    if (std::find(branchNames.begin(), branchNames.end(), "KE") != branchNames.end())
        tree->SetBranchAddress("KE", &dataKE);
    else
        Fatal("collectData", "Requested branch 'E' does not exist in the file.");

    // Get the number of entries
    int entries = tree->GetEntries();

    // Collect the data
    for (int i = 0; i < entries; i++) {
        // Get the corresponding entry
        tree->GetEntry(i);

        // Collect the data
        if (virtualdetectorId == dataVirtualdetectorId) {
            pdgIds.push_back(dataPdgId);
            energies.push_back(dataKE);
        };
    };

    // Check if data has been collected
    if (energies.size() == 0)
        Fatal("collectData", "No data was collected from this file");

    double eMin = *std::min_element(energies.begin(), energies.end());
    if (eMin < (-1 * std::numeric_limits<double>::epsilon()))
        Fatal("collectVirtualDetectorData", "Minimum of the energy is negative, there is an issue with the table generation code");
    // Clear
    file->Close();
    delete file;
    std::cout << "Finished processing file " << fileName << std::endl;
    return;
};

void printTable(std::vector<int> &pdgIds, const std::vector<double> &energies, const long long nPOTs, const double &virtualdetectorRadius, const long long &numBatchesPerSuperCycle, const int w = 95) {
    /*
        Description
        Calculates the parameter values for the number of particles, average particle energy, energy flux, and POT normalized intensity

        Arguments
        pdgIds - vector of relevant PDG IDs
        energies - vector of relevant particle energies, in MeV
        nPOTs - number of protons of target equivalent for the given dataset
        virtualdetectorRadius - as documented in function "table"
        numBatchesPerSuperCycle - as documented in function "table"
        w - as documented in function "table"

        Variables
        pdgIdsSet - collection of unique PDG IDs
        nPdgIds - number of unique PDG IDs
        resultCount - vector of particle count by PDG ID
        resultE - vector of total energy by PDG ID
        nEnergies - number of entries in the tree
        pdgId - entry PDG ID
        index - position of PDG ID in pdgIdsSet
        highEnergyPhotonCount - count of high energy photons
        highEnergyPhotonE - total energy of high energy photons
        highEnergyPhotonMinE - energy threshold of high energy photons
        it - iterator for pdgIdSet
        nPOTsPerMicroSpill - number of POTs per micro spill
        nMicroSpillsPerMacroSpill - number of micro spills per macro spill
        nMacroSpillsPerSuperCycle - number of macro spills per super cycle
        nMicroSpills - total number of complete micro spills based on the POT count
        nMacroSpills - total number of complete macro spills based on the micro spill count
        nSuperCycles - total number of complete super cycles spills based on the macro spill count
        nRemPOTs - number of POTs in incomplete micro spills
        nRemMicroSpills - number of micro spills in incomplete macro spills
        nRemMacroSpills - number of macro spills in incomplete super cycles
        tMicroSpill - duration of the micro spill, in seconds
        tMacroSpill - duration of the macro spill, in seconds
        tMacroSpill45 - duration of the extra break betweeen spill 4 and spill 5, in seconds
        tSuperCycle - duration of the super cycle, in seconds
        time - equivalent time of beam operations based on the POT count and the booster batch count, in seconds
        i - iterator for pdgIdsSet
        stream - string stream for converting max photon energy to a string
        highEnergyPhotonTitle - high energy photon title for the table
    */
    // Set up data collation variables
    std::set<int> pdgIdsSet(pdgIds.begin(), pdgIds.end());
    int nPdgIds = pdgIdsSet.size();
    std::vector<int>            resultCount(nPdgIds, 0);
    std::vector<long double>    resultE(nPdgIds, 0);

    // Collate the data
    int nEnergies = energies.size(), pdgId = 0, index = 0, highEnergyPhotonCount = 0;
    double highEnergyPhotonE = 0, E = 0;
    const double highEnergyPhotonMinE = 0.1;
    for (int i = 0; i < nEnergies; i++) {
        pdgId = pdgIds[i];
        E = energies[i];
        index = std::distance(pdgIdsSet.begin(), pdgIdsSet.find(pdgId));
        resultCount[index]++;
        resultE[index] += E;
        if (pdgId == 22 && E > highEnergyPhotonMinE) {
            highEnergyPhotonCount++;
            highEnergyPhotonE += E;
        }
    };

    // Determine the effective real time from nPOTs
    std::set<int>::iterator it = pdgIdsSet.begin();
    const int nPOTsPerMicroSpill         = (numBatchesPerSuperCycle == 2 ? 39e6 : 2e7),     nMicroSpills = nPOTs        / nPOTsPerMicroSpill,           nRemPOTs        = nPOTs        % nPOTsPerMicroSpill;
    const int nMicroSpillsPerMacroSpill  = 31858,                                           nMacroSpills = nMicroSpills / nMicroSpillsPerMacroSpill,    nRemMicroSpills = nMicroSpills % nMicroSpillsPerMacroSpill;
    const int nMacroSpillsPerSuperCycle  = (numBatchesPerSuperCycle == 2 ? 8 : 4),          nSuperCycles = nMacroSpills / nMacroSpillsPerSuperCycle,    nRemMacroSpills = nMacroSpills % nMacroSpillsPerSuperCycle;
    const long double tMicroSpill = 1695e-9, tMacroSpill = 59e-3, tMacroSpill45 = 31e-3, tSuperCycle = 1.33;

    const long double time =  nSuperCycles * tSuperCycle                                                // Count the number of supercycle intervals
                            + nRemMacroSpills * tMacroSpill                                             // Add the remaining macrospills
                            + (nMacroSpills > 3 && numBatchesPerSuperCycle == 2 ? tMacroSpill45 : 0)    // Add the interval between the 4th and 5th macrospills, if relevant
                            + nRemMicroSpills * tMicroSpill                                             // Add the remaining microspills
                            + tMicroSpill * (nRemPOTs / nPOTsPerMicroSpill);                            // Add the remaining POTs

    // Determine the virtual detector area
    const long double pi = 3.141592653589793, virtualdetectorArea = pi * virtualdetectorRadius * virtualdetectorRadius;

    // Set up the table
    std::cout << std::string(w, '-')   << std::endl; // Title line
    std::cout << std::setw(15)          << std::left << "PDG ID"
              << std::setw(10)          << "N"
              << std::setw(15)          << "Avg E [MeV]"
              << std::setw(25)          << "E Flux [TeV/cm2/s]"
              << std::setw(30)          << "POT norm. intensity [/cm2/POT]" << std::endl;
    std::cout << std::string(w, '-')   << std::endl; // Separator line

    // Print the table
    for (int i = 0; i < nPdgIds; i++) {
        std::cout << std::setw(15) << std::scientific << std::setprecision(2) << std::left << *it
                  << std::setw(10) << std::scientific << std::setprecision(2) << resultCount[i]
                  << std::setw(15) << std::scientific << std::setprecision(2) << resultE[i]/resultCount[i]
                  << std::setw(25) << std::scientific << std::setprecision(2) << resultE[i]/(1e6 * virtualdetectorArea * time)
                  << std::setw(30) << std::scientific << std::setprecision(2) << resultCount[i]/(virtualdetectorArea * nPOTs) << std::endl;
        std::advance(it, 1);
    };
    // Print the high energy photon information
    std::stringstream stream;
    stream << std::fixed << std::setprecision(2) << highEnergyPhotonMinE;
    const std::string highEnergyPhotonTitle = std::string("E > ") + stream.str() + " MeV";
    std::prev(it);
    std::cout << std::setw(15) << std::left << highEnergyPhotonTitle
              << std::setw(10) << highEnergyPhotonCount
              << std::setw(15) << highEnergyPhotonE/highEnergyPhotonCount
              << std::setw(25) << highEnergyPhotonE/(1e6 * virtualdetectorArea * time) // 1e6 converts MeV to TeV
              << std::setw(30) << highEnergyPhotonE/(virtualdetectorArea * nPOTs) << std::endl;
    std::cout << std::string(w, '-')   << std::endl; // Finishing line
    std::cout << std::endl; // Clearance line

    return;
};

void table(const std::vector<std::string> electronFileNames, const long long electronPOTs, const std::vector<std::string> muonFileNames, const long long muonPOTs, const std::string treeName, const int virtualdetectorId, const double virtualdetectorRadius, const int numBatchesPerSuperCycle) {
    /*
        Description
        Collate the data determined in STM studies to generate summary table
        Uses the results of MakeTree.fcl
        usage example - see table.sh

        Arguments
        electronFileNames - vector of EleBeamCat derived ROOT file names as a relative path to cwd
        electronPOTs - number of POTs used to generate all the files listed in electronFileNames
        muonFileNames - vector of MuBeamCat derived ROOT file names as a relative path to cwd
        muonPOTs - number of POTs used to generate all the files listed in muonFileNames
        treeName - name of tree containing all the data for electronFileNames and muonFileNames
        virtualdetectorId - Virtual Detector ID for plotting
        virtualdetectorRadius - radius of selected Virtual Detector, in cm
        numBatchesPerSuperCycle - number of batches per supercycle either 1 for single batch mode or 2 for dual batch mode

        Variables
        allowedBatchModes - used to check the range of numBatchesPerSuperCycle, in MeV
        electronEnergies - vector for energies from EleBeamCat results
        muonEnergies - vector of energies from MuBeamCat results, in MeV
        electronPdgIds - vector of PDG IDs from EleBeamCat results
        muonPdgIds - vector of PDG IDs from MuBeamCat results
        fileName - iterator through electronFileNames and muonFileNames
        w - width of the table in chars
        spacing - used to format the table header
    */
    SetErrorHandler(customErrorHandler);

    // Pre data handling checks
    if (!electronPOTs || !muonPOTs)
        Fatal("table", "Number of POTs associated with both the EleBeamCat and MuBeamCat results must be non zero");
    if (treeName == "")
        Fatal("table", "No tree name provided");
    if (!virtualdetectorId)
        Fatal("table", "Virtual detector ID cannot be 0");
    if (virtualdetectorRadius < std::numeric_limits<double>::epsilon())
        Fatal("table", "Virtual detector radius cannot be 0");
    std::vector<int> allowedBatchModes = {1, 2};
    if (std::find(allowedBatchModes.begin(), allowedBatchModes.end(), numBatchesPerSuperCycle) == allowedBatchModes.end())
        Fatal("table", "Number of batches per super cycle has to be either 1 or 2");

    // Set up data collection variables
    std::vector<double> electronEnergies,   muonEnergies;
    std::vector<int>    electronPdgIds,     muonPdgIds;

    // Collect data
    for (const std::string fileName : electronFileNames)
        collectVirtualdetectorData(fileName, treeName, virtualdetectorId, electronEnergies, electronPdgIds);
    for (const std::string fileName : muonFileNames)
        collectVirtualdetectorData(fileName, treeName, virtualdetectorId, muonEnergies,     muonPdgIds);

    // Validate that data was found
    if (electronEnergies.empty() && muonEnergies.empty())
        Fatal("table", "No data found");
    std::cout << "\n" << std::endl; // Clear line to sepaate data collection and summarizing

    int w = 95, spacing = 0;
    std::string header = "";
    // Generate EleBeamCat table
    if (!electronEnergies.empty()) {
        header = "EleBeamCat at VD" + std::to_string(virtualdetectorId);
        spacing = (w - header.length()) / 2;
        std::cout << std::string(spacing, ' ') << header << std::endl;
        printTable(electronPdgIds,  electronEnergies,   electronPOTs,   virtualdetectorRadius, numBatchesPerSuperCycle);
    }
    else
        std::cout << "Not generating EleBeamCat table, no data found" << std::endl;

    // Generate MuBeamCat table
    if (!muonEnergies.empty()) {
        header = "MuBeamCat at VD" + std::to_string(virtualdetectorId);
        spacing = (w - header.length()) / 2;
        std::cout << std::string(spacing, ' ') << header << std::endl;
        printTable(muonPdgIds,      muonEnergies,       muonPOTs,        virtualdetectorRadius, numBatchesPerSuperCycle);
    }
    else
        std::cout << "Not generating MuBeamCat table, no data found" << std::endl;
    return;
};
