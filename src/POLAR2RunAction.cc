#include "POLAR2RunAction.hh"
#include "POLAR2Run.hh"
#include "POLAR2GlobalConfig.hh"
#include "G4RunManager.hh"
#include "Randomize.hh"
#include "g4root.hh"
#include "Randomize.hh"
#include <sstream>

POLAR2RunAction::POLAR2RunAction(G4String the_output_file, G4bool fixed_name) {
    output_file_ = the_output_file;
    fixed_output_ = fixed_name;
}

POLAR2RunAction::~POLAR2RunAction() {

}

G4Run* POLAR2RunAction::GenerateRun() {
    return new POLAR2Run();
}

G4int POLAR2RunAction::randomSeed_ = 0;
G4int POLAR2RunAction::runId_ = 0;

void POLAR2RunAction::BeginOfRunAction(const G4Run* aRun) {
    G4AnalysisManager* fAnalysisManager = G4AnalysisManager::Instance();

    if (IsMaster()) {
        randomSeed_ = G4Random::getTheSeed();
        runId_      = aRun->GetRunID();
    }
    G4String filename(output_file_);
    // append extra info to filename
    if (!fixed_output_) {
        G4String extension;
        if (filename.find(".") != G4String::npos) {
            extension = filename.substr(filename.find("."));
            filename = filename.substr(0, filename.find("."));
        } else {
            extension = "";
        }
        std::ostringstream os;
        os << "_s" << randomSeed_ << "_r" << runId_;
        filename.append(os.str());
        filename.append(extension);
    }
    // open file
    POLAR2GlobalConfig* fPOLAR2GlobalConfig = POLAR2GlobalConfig::Instance();
    if (fPOLAR2GlobalConfig->output_directory.size() > 0) {
        filename.prepend("/");
        filename.prepend(fPOLAR2GlobalConfig->output_directory.c_str());
    }
    fAnalysisManager->OpenFile(filename);

}

void POLAR2RunAction::EndOfRunAction(const G4Run* aRun) {
    G4AnalysisManager* fAnalysisManager = G4AnalysisManager::Instance();

    const POLAR2Run* theRun = static_cast<const POLAR2Run*>(aRun);
    if (IsMaster()) {
        // fill EndOfRun information
        fAnalysisManager->FillNtupleIColumn(0, 0, theRun->GetRunID());
        fAnalysisManager->FillNtupleIColumn(0, 1, theRun->GetNumberOfEvent());
        fAnalysisManager->FillNtupleIColumn(0, 2, theRun->GetTotalHits());
        fAnalysisManager->FillNtupleIColumn(0, 3, theRun->GetTotalEvents());
        fAnalysisManager->AddNtupleRow(0);
    }

    fAnalysisManager->Write();
    fAnalysisManager->CloseFile();
}
