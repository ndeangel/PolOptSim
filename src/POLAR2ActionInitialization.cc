#include "POLAR2ActionInitialization.hh"
#include "POLAR2PrimaryGeneratorAction.hh"
#include "POLAR2RunAction.hh"
#include "POLAR2EventAction.hh"
#include "POLAR2SteppingAction.hh"
#include "g4root.hh"
#include "G4RunManager.hh"
#include "POLAR2GlobalConfig.hh"

POLAR2ActionInitialization::POLAR2ActionInitialization(G4bool the_gps_flag, G4String the_output_file, G4bool fixed_name) {
    fParticleDataFileR_  = new ParticleDataFileR();
    gps_flag_    = the_gps_flag;
    output_file_ = the_output_file;
    fixed_filename_ = fixed_name;
}

POLAR2ActionInitialization::~POLAR2ActionInitialization() {
    delete fParticleDataFileR_;
}

ParticleDataFileR* POLAR2ActionInitialization::GetParticleDataFileR() {
    return fParticleDataFileR_;
}

void POLAR2ActionInitialization::Build() const {
    SetUserAction(new POLAR2PrimaryGeneratorAction(gps_flag_, fParticleDataFileR_));
    SetUserAction(new POLAR2RunAction(output_file_, fixed_filename_));
    SetUserAction(new POLAR2EventAction());
    SetUserAction(new POLAR2SteppingAction());

    POLAR2GlobalConfig* fPOLAR2GlobalConfig = POLAR2GlobalConfig::Instance();

    // prepare file
    G4AnalysisManager* fAnalysisManager = G4AnalysisManager::Instance();
    if (G4RunManager::GetRunManager()->GetRunManagerType() == G4RunManager::RMType::sequentialRM) {
        fAnalysisManager->CreateNtuple("Keyword", "Keyword");
        fAnalysisManager->CreateNtupleIColumn("RunID");                     // 0
        fAnalysisManager->CreateNtupleIColumn("NumberOfEvents");            // 1
        fAnalysisManager->CreateNtupleIColumn("TotalHits");                 // 2
        fAnalysisManager->CreateNtupleIColumn("TotalEvents");               // 3
        fAnalysisManager->FinishNtuple();
    }

    fAnalysisManager->CreateNtuple("Primary", "Primary");
    // necessary
    fAnalysisManager->CreateNtupleIColumn("EventID");                       // 0
    fAnalysisManager->CreateNtupleIColumn("HitsCount");                     // 1
    // optional
    if (fPOLAR2GlobalConfig->simout_more) {
        fAnalysisManager->CreateNtupleSColumn("PName_Primary");             // 2
        fAnalysisManager->CreateNtupleFColumn("X_Primary");                 // 3
        fAnalysisManager->CreateNtupleFColumn("Y_Primary");                 // 4
        fAnalysisManager->CreateNtupleFColumn("Z_Primary");                 // 5
        fAnalysisManager->CreateNtupleDColumn("Time_Primary");              // 6
        fAnalysisManager->CreateNtupleFColumn("Energy_Primary");            // 7
        fAnalysisManager->CreateNtupleFColumn("POLX_Primary");              // 8
        fAnalysisManager->CreateNtupleFColumn("POLY_Primary");              // 9
        fAnalysisManager->CreateNtupleFColumn("POLZ_Primary");              // 10
        fAnalysisManager->CreateNtupleFColumn("Theta_Primary");             // 11
        fAnalysisManager->CreateNtupleFColumn("Phi_Primary");               // 12
    }
    fAnalysisManager->FinishNtuple();

    fAnalysisManager->CreateNtuple("HitsCol", "HitsCollection");
    // necessary
    fAnalysisManager->CreateNtupleIColumn("EventID");                       // 0
    fAnalysisManager->CreateNtupleIColumn("TrackID");                       // 1
    fAnalysisManager->CreateNtupleIColumn("ParticleCode");                  // 2
    fAnalysisManager->CreateNtupleIColumn("BarID");                         // 3
    fAnalysisManager->CreateNtupleIColumn("ModID");                         // 4
    fAnalysisManager->CreateNtupleFColumn("LocPosX");                       // 5
    fAnalysisManager->CreateNtupleFColumn("LocPosY");                       // 6
    fAnalysisManager->CreateNtupleFColumn("LocPosZ");                       // 7
    fAnalysisManager->CreateNtupleFColumn("EnergyVis");                     // 8
    // optional
    if (fPOLAR2GlobalConfig->simout_more) {
        fAnalysisManager->CreateNtupleSColumn("ParticleName");              // 9
        fAnalysisManager->CreateNtupleDColumn("PreGlobalTime");             // 10
        fAnalysisManager->CreateNtupleIColumn("StepStatus");                // 11
        fAnalysisManager->CreateNtupleFColumn("PreMomTheta");               // 12
        fAnalysisManager->CreateNtupleFColumn("PreMomPhi");                 // 13
        fAnalysisManager->CreateNtupleFColumn("PostMomTheta");              // 14
        fAnalysisManager->CreateNtupleFColumn("PostMomPhi");                // 15
        fAnalysisManager->CreateNtupleFColumn("PreKinEnergy");              // 16
        fAnalysisManager->CreateNtupleFColumn("PostKinEnergy");             // 17
        fAnalysisManager->CreateNtupleSColumn("ProcessName");               // 18
        fAnalysisManager->CreateNtupleFColumn("EnergyDep");                 // 19
        fAnalysisManager->CreateNtupleDColumn("DeltaTime");                 // 20
    }
    fAnalysisManager->FinishNtuple();

}

void POLAR2ActionInitialization::BuildForMaster() const {
    SetUserAction(new POLAR2RunAction(output_file_, fixed_filename_));
    // prepare file
    G4AnalysisManager* fAnalysisManager = G4AnalysisManager::Instance();
    fAnalysisManager->CreateNtuple("Keyword", "Keyword");
    fAnalysisManager->CreateNtupleIColumn("RunID");
    fAnalysisManager->CreateNtupleIColumn("NumberOfEvents");
    fAnalysisManager->CreateNtupleIColumn("TotalHits");
    fAnalysisManager->CreateNtupleIColumn("TotalEvents");
    fAnalysisManager->FinishNtuple();

}
