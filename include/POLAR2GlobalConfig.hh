#ifndef POLAR2GlobalConfig_HH
#define POLAR2GlobalConfig_HH

#include "globals.hh"

class POLAR2ConfigMessenger;

class POLAR2GlobalConfig {
private:
    POLAR2ConfigMessenger* fPOLAR2ConfigMessenger_;

public:
    POLAR2GlobalConfig();
    ~POLAR2GlobalConfig();

public:
    static POLAR2GlobalConfig* Instance();

private:
    static POLAR2GlobalConfig* fgInstance_;

public:
    // config parameter list
    G4double hit_threshold;  // lowest energy with unit
    G4double bar_threshold;
    G4String output_directory;
    G4double polarisation_angle;
    G4double incidence_theta;
    G4double incidence_phi;
    G4int    event_verbose;
    G4bool   primary_only;
    G4double birks_constant;
    G4bool   spacelab;
    G4int    phys_verbose;
    G4bool   simout_more;
    G4bool   phys_more;
    G4double antenna_angle_ud;
    G4double antenna_angle_lr;

public:
    void print_config();

};

#endif
