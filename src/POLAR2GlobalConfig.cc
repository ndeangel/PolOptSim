#include "POLAR2GlobalConfig.hh"
#include "POLAR2ConfigMessenger.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

POLAR2GlobalConfig* POLAR2GlobalConfig::fgInstance_ = NULL;

POLAR2GlobalConfig::POLAR2GlobalConfig() {
    if (fgInstance_ != NULL) {
        G4ExceptionDescription description;
        description << "      "
                    << "POLAR2GlobalConfig already exists."
                    <<  "Cannot create another instance.";
        G4Exception("POLAR2GlobalConfig::POLAR2GlobalConfig()", "Config_F001",
                FatalException, description);
    }
    fPOLAR2ConfigMessenger_ = new POLAR2ConfigMessenger(this);
    fgInstance_ = this;

    // initial value
    hit_threshold = 0.0001*eV;
    bar_threshold = 0.0001*eV;
    output_directory = "";
    polarisation_angle = 0.0*deg;
    incidence_theta = 0.0*deg;
    incidence_phi = 0.0*deg;
    event_verbose = 0;
    primary_only = false;
    birks_constant = 0.143;  //0.143
    spacelab = false;
    phys_verbose = 0;
    simout_more = false;
    phys_more = false;
    antenna_angle_ud = 0.0*deg;
    antenna_angle_lr = 0.0*deg;
}

POLAR2GlobalConfig* POLAR2GlobalConfig::Instance() {
    if (fgInstance_ == NULL) {
        fgInstance_ = new POLAR2GlobalConfig();
    }
    return fgInstance_;
}

POLAR2GlobalConfig::~POLAR2GlobalConfig() {
    delete fPOLAR2ConfigMessenger_;
}

void POLAR2GlobalConfig::print_config() {
    G4cout << "=============================================================" << G4endl;
    G4cout << "======================= Configuration =======================" << G4endl;
    G4cout << "=============================================================" << G4endl;
    G4cout << " - hit_threshold        = " << hit_threshold / eV << " eV" << G4endl;
    G4cout << " - bar_threshold        = " << bar_threshold / eV << " eV" << G4endl;
    G4cout << " - output_directory     = " << output_directory << G4endl;
    G4cout << " - polarisation_angle   = " << polarisation_angle / deg << " deg" << G4endl;
    G4cout << " - incidence_theta      = " << incidence_theta / deg << " deg" << G4endl;
    G4cout << " - incidence_phi        = " << incidence_phi / deg << " deg" << G4endl;
    G4cout << " - event_verbose        = " << event_verbose << G4endl;
    G4cout << " - primary_only         = " << primary_only << G4endl;
    G4cout << " - birks_constant       = " << birks_constant << " mm/MeV" << G4endl;
    G4cout << " - spacelab             = " << spacelab << G4endl;
    G4cout << " - phys_verbose         = " << phys_verbose << G4endl;
    G4cout << " - simout_more          = " << simout_more << G4endl;
    G4cout << " - phys_more            = " << phys_more << G4endl;
    G4cout << " - antenna_angle_ud     = " << antenna_angle_ud / deg << " deg" << G4endl;
    G4cout << " - antenna_angle_lr     = " << antenna_angle_lr / deg << " deg" << G4endl;
    G4cout << "=============================================================" << G4endl;
}
