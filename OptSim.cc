#include "G4RunManager.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#endif

#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#include "G4VisExecutive.hh"

#include "POLAR2DetectorConstruction.hh"
#include "POLAR2PhysicsList.hh"
#include "POLAR2ActionInitialization.hh"
#include "POLAR2GlobalConfig.hh"

#include "OptionsManager.hh"

int main(int argc,char** argv) {
    // process command parameters
    OptionsManager options_mgr;
    if (!options_mgr.parse(argc, argv)) {
        if (options_mgr.get_version_flag()) {
            options_mgr.print_version();
        } else {
            options_mgr.print_help();
        }
        return 2;
    }

    POLAR2GlobalConfig* fPOLAR2GlobalConfig = POLAR2GlobalConfig::Instance();

    G4UImanager* uiManager = G4UImanager::GetUIpointer();
    G4String execute = "/control/execute ";
    if (options_mgr.config_file != "") {
        uiManager->ApplyCommand(execute + options_mgr.config_file);
        fPOLAR2GlobalConfig->print_config();
    }

    POLAR2DetectorConstruction* fPOLAR2DetectorConstruction = new POLAR2DetectorConstruction();
//     POLAR2PhysicsList*          fPOLAR2PhysicsList          = new POLAR2PhysicsList(fPOLAR2GlobalConfig->phys_verbose);
    POLAR2ActionInitialization* fPOLAR2ActionInitialization = new POLAR2ActionInitialization(options_mgr.gps_flag, options_mgr.output_file, options_mgr.fixed_name);

#ifdef G4MULTITHREADED
    G4MTRunManager* runManagerMT = NULL;
    G4RunManager*   runManagerSQ = NULL;
    if (options_mgr.mt_flag && options_mgr.gps_flag) {
        runManagerMT = new G4MTRunManager();
        G4cout << "#### using multi-threaded mode with " << options_mgr.num_of_thread << " threads ####" << G4endl;
        runManagerMT->SetNumberOfThreads(options_mgr.num_of_thread);
        // set mandatory initialization classes
        runManagerMT->SetUserInitialization(fPOLAR2DetectorConstruction);
//         runManagerMT->SetUserInitialization(fPOLAR2PhysicsList);
	runManagerMT->SetUserInitialization(new PhysicsList);

        runManagerMT->SetUserInitialization(fPOLAR2ActionInitialization);
        // initialize G4 kernel
        G4cout << "initialize G4 kernel" << G4endl;
        runManagerMT->Initialize();
    } else {
        runManagerSQ = new G4RunManager();
        G4cout << "#### using sequential mode ####" << G4endl;
        // set mandatory initialization classes
        runManagerSQ->SetUserInitialization(fPOLAR2DetectorConstruction);
        runManagerSQ->SetUserInitialization(fPOLAR2PhysicsList);
        runManagerSQ->SetUserInitialization(fPOLAR2ActionInitialization);
        // initialize G4 kernel
        G4cout << "initialize G4 kernel" << G4endl;
        runManagerSQ->Initialize();
    }
#else
    G4RunManager* runManager = new G4RunManager();
    G4cout << "#### WARNING: only sequential mode can be used ####" << G4endl;
    G4cout << "#### using sequential mode ####" << G4endl;
    // set mandatory initialization classes
    runManager->SetUserInitialization(fPOLAR2DetectorConstruction);
//     runManager->SetUserInitialization(fPOLAR2PhysicsList);
    runManager->SetUserInitialization(new PhysicsList);

    runManager->SetUserInitialization(fPOLAR2ActionInitialization);
    // initialize G4 kernel
    G4cout << "initialize G4 kernel" << G4endl;
    runManager->Initialize();
#endif

    G4cout << "#### G4 KERNEL INITIALIZED ####" << G4endl;

    G4VisManager* visManager = NULL;
    if (options_mgr.gui_flag || options_mgr.ter_flag) {
        visManager = new G4VisExecutive("warnings"); // "quiet", "errors", "warnings"
        visManager->Initialize();
        G4cout << "#### G4 VIS INITIALIZED ####" << G4endl;
    }


    if (options_mgr.gui_flag) {
        G4UIExecutive* ui = new G4UIExecutive(1, argv);
        uiManager->ApplyCommand(execute + options_mgr.vis_mac_file);
        ui->SessionStart();
        delete ui;

    } else if (options_mgr.ter_flag) {
        G4UIsession* session = NULL;
#ifdef G4UI_USE_TCSH
        session = new G4UIterminal(new G4UItcsh);
#else
        session = new G4UIterminal();
#endif
        uiManager->ApplyCommand(execute + options_mgr.vis_mac_file);
        session->SessionStart();
        delete session;
    } else {
        // start simulation
        G4cout<< "************** good morning **************" << G4endl;

        if (options_mgr.gps_flag) {
            G4cout << "==== using general particle source ====" << G4endl;
            uiManager->ApplyCommand(execute + options_mgr.primary_file);
        } else {
            G4cout << "==== using root file of particle data ====" << G4endl;
            if (!fPOLAR2ActionInitialization->GetParticleDataFileR()->open(options_mgr.primary_file.c_str())) {
                G4cout << "particle data file open failed." << G4endl;
                return 1;
            }
        }

        uiManager->ApplyCommand(execute + options_mgr.run_mac_file);

        if (!options_mgr.gps_flag) {
            fPOLAR2ActionInitialization->GetParticleDataFileR()->close();
        }

        G4cout<< "************** good night **************" << G4endl;
    }

    // job termination
    if (options_mgr.gui_flag || options_mgr.ter_flag) {
        delete visManager;
    }
#ifdef G4MULTITHREADED
    if (options_mgr.mt_flag) {
        delete runManagerMT;
    } else {
        delete runManagerSQ;
    }
#else
    delete runManager;
#endif
    return 0;
}

