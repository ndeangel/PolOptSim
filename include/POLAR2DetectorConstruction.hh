#ifndef POLAR2DetectorConstruction_h
#define POLAR2DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"

//Unit
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4Material.hh"
#include "G4OpticalSurface.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"

class G4Material;

class POLAR2DetectorConstruction : public G4VUserDetectorConstruction {
public:
POLAR2DetectorConstruction();
~POLAR2DetectorConstruction();

public:
G4VPhysicalVolume* Construct();
void ConstructSDandField();

private: // Materials
bool materials_defined_;
G4Material* Environment_;
G4Material* theScint;
G4Material* CarbonFiber_;
G4Material* Steel_;
G4Material* detWMaterial;
G4Material* GasMaterial;
G4Material* Rubber_;
//G4Material* claryl;
//G4Material* reflcoating;
G4Material* PEEK_Plastic_;
G4Material* G4_Mag_;

G4Material* G4_Al_;
G4Material* G4_Cu_;
G4Material* PMTdynodeMat_;
G4Material* VikuitiESR_;
G4Material* Polythene_Al_;
G4Material* Resin_;
G4Material* FR4_;
G4Material* Sylgard_184_;
G4Material* AlloyAl_7075_;
G4Material* AlloyAl_2219_;
G4Material* Accura25;
G4Material* dc93_500;
G4Material* EpoxyResin;

G4OpticalSurface* OpSiPMResinSurface;
G4OpticalSurface* OpResinSiPMSurface;

G4OpticalSurface* OpPadResinSurface;
G4OpticalSurface* OpResinPadSurface;

G4OpticalSurface* OpAirWrapSurface;
G4OpticalSurface* OpWrapAirSurface;
G4OpticalSurface* OpAirGridSurface;
G4OpticalSurface* OpGridAirSurface;
G4OpticalSurface* OpSciAirSurface;
G4OpticalSurface* OpAirSciSurface;

G4OpticalSurface* OpAirPerfReflSurface;
G4OpticalSurface* OpPerfReflAirSurface;

G4OpticalSurface* OpPadSciSurface;
G4OpticalSurface* OpSciPadSurface;

G4OpticalSurface* OpAirTeflonSurface; //LEAP
G4OpticalSurface* OpTeflonAirSurface;


G4int nScintBars = 64;
G4double BarPitch = 6.2*mm;

G4double Scint_length = 125.0*mm;
G4double Scint_thick = 5.9*mm;
G4double Grid_height = 3.0*mm;
G4double Vikuiti_thickness = 0.065*mm;

G4double NoWrap_height = 3.5*mm;
G4double Coating_height = 4.0*mm;

//G4double Scint_end = 5.0*mm;
//G4double Trunc_height= 5.0*mm;

G4double Pad_height = 0.3*mm; //0.3*mm
G4double EpoxyResin_height = 0.1*mm;
G4double SiPMsilicon_height = 1.2*mm;

//bool use_coating = true;
bool use_EJ200 = false;
bool use_SiPM_PDE = false;  // whether to use the SiPM PDE or a perfect 100% absorbing detector
bool vikuiti_is_metal = true; // whether to use dielectric_metal or dielectric_dielectric for the Vikuiti-Air optical surfaces

G4double sigma_alpha = 1.82*deg;  //3.45  1.82

G4double alutick = 0.01*mm;
G4double coat_thick = 0.005*mm;
G4double Air_thick_wrap = 0.01*mm;
G4double Air_thick_coat = 0.045*mm;


// LEAP
// whether to use a conical tapering for the LEAP scintillator
bool LEAPScint_conical = true;

bool LEAPvisualization = false;
G4double LEAPScint_length = 100.0*mm;
G4double LEAPScint_thick = 17.0*mm;  // 17.0*mm, tried with 8.5*mm "small square"
G4double LEAP_PMTEntranceRadius = 4.0*mm;
G4double LEAP_PMTRadius = 7.0*mm;
G4double LEAP_PMT_thick = 10.0*mm;
G4double LEAP_Optpad_apperture_thick = 0.75*mm;
G4double LEAP_Optpad_thick = 0.8*mm;
G4double LEAPTefloncone_thick = 12*mm;
G4double LEAPTeflonsquare_thick = 4.5*mm;


private: // Volumes
G4LogicalVolume* WorldLog_;

G4LogicalVolume* ScintLog_;
G4LogicalVolume* SiPMLog_;

G4LogicalVolume* PerfRefl_x_log;
G4LogicalVolume* PerfRefl_y_log;
G4LogicalVolume* PerfRefl_z_log;

G4LogicalVolume* vikuit_x_log;
G4LogicalVolume* vikuit_y_log;
G4LogicalVolume* vikuit_z_log;

G4LogicalVolume* air_inner_x_log;
G4LogicalVolume* air_outer_x_log;
G4LogicalVolume* air_inner_y_log;
G4LogicalVolume* air_outer_y_log;
G4LogicalVolume* air_z_log;

G4LogicalVolume* GridLog_;
G4LogicalVolume* l_grid_long_thin;
G4LogicalVolume* l_grid_long_thick;
G4LogicalVolume* l_grid_short_thin;
G4LogicalVolume* l_grid_short_thick;
G4LogicalVolume* l_cross_long;
G4LogicalVolume* l_cross_short;

G4LogicalVolume* opt_pad_log;
G4LogicalVolume* EpoxyResin_log;

G4LogicalVolume* reflcoating_x_log;
G4LogicalVolume* reflcoating_y_log;

G4LogicalVolume* claryl_x_log;
G4LogicalVolume* claryl_y_log;
G4LogicalVolume* claryl_z_log;

G4LogicalVolume* DetectorLog_;

G4MaterialPropertiesTable* SciMPT;

G4LogicalVolume* air_test_log;
G4LogicalVolume* airbbox;
G4LogicalVolume* grid_plate_test;
G4LogicalVolume* ScintBlockLog_;

//LEAP
G4LogicalVolume* LEAPScintLog_;
G4LogicalVolume* LEAPCylindricalScintLog_;
G4LogicalVolume* LEAPvikuit_x_log;
G4LogicalVolume* LEAPvikuit_y_log;
G4LogicalVolume* LEAPvikuit_z_log;
G4LogicalVolume* LEAPCylindricalvikuit_xy_log;
G4LogicalVolume* LEAPCylindricalvikuit_z_log;
G4LogicalVolume* LEAPair_x_log;
G4LogicalVolume* LEAPair_y_log;
G4LogicalVolume* LEAPair_z_log;
G4LogicalVolume* LEAPCylindricalair_xy_log;
G4LogicalVolume* LEAPCylindricalair_z_log;
G4LogicalVolume* LEAPopt_pad_log;
G4LogicalVolume* LEAPPMTLog_;
G4LogicalVolume* TeflonLog_;
G4LogicalVolume* LEAPair_teflonsquare_log;
G4LogicalVolume* LEAPair_teflonconical_log;

private:
void DefineMaterials_();
void ConstructModule_();
void ConstructDetector_();
void ConstructSpaceLab_();


};

#endif /* POLAR2DetectorConstruction_h */
