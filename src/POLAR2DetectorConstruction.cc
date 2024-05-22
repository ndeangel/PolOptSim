#include "POLAR2DetectorConstruction.hh"

//Unit
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

//material
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4NistManager.hh"

//geometry
#include "G4Orb.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Trap.hh"
#include "G4Polycone.hh"
#include "G4Polyhedra.hh"
#include "G4Sphere.hh"
#include "G4Para.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4MultiUnion.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include <cstdio>
#include <cmath>
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4Trd.hh"

//transformation
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"

//sensitive detector
#include "G4SDManager.hh"
#include "POLAR2SensitiveDetector.hh"

//visual attribute
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

//global function
#include "globals.hh"

// configure
#include "POLAR2GlobalConfig.hh"

POLAR2DetectorConstruction::POLAR2DetectorConstruction() {
    WorldLog_          = NULL;
    ScintLog_		= NULL;
    materials_defined_ = false;
}

POLAR2DetectorConstruction::~POLAR2DetectorConstruction() {

}

G4VPhysicalVolume* POLAR2DetectorConstruction::Construct() {
    //-----------------------------------------------------------------------------------------------------------
    //--------Define the detector geometry
    //-----------------------------------------------------------------------------------------------------------

    POLAR2GlobalConfig* fPOLAR2GlobalConfig = POLAR2GlobalConfig::Instance();

    DefineMaterials_();

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    // Construct WorldPhys
    G4double world_hx = (fPOLAR2GlobalConfig->spacelab ? 7.5*m : 500*mm);
    G4double world_hy = (fPOLAR2GlobalConfig->spacelab ? 7.5*m : 500*mm);
    G4double world_hz = (fPOLAR2GlobalConfig->spacelab ? 7.5*m : 500*mm);
    G4Box* WorldBox = new G4Box("WorldBox", world_hx, world_hy, world_hz);
    G4LogicalVolume* WorldLog = new G4LogicalVolume(WorldBox, Environment_, "WorldLog");
    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    if (!fPOLAR2GlobalConfig->primary_only) {
        // Construct Detector
        ConstructModule_();
        ConstructDetector_();
        new G4PVPlacement(NULL, G4ThreeVector(), DetectorLog_, "Detector", WorldLog, false, 0);
        // Construct SpaceLab
//         if (fPOLAR2GlobalConfig->spacelab) {
//             G4cout << "Building SpaceLab TG02 ..." << G4endl;
//             ConstructSpaceLab_();
//             new G4PVPlacement(NULL, G4ThreeVector(), SpaceLabLog_, "SpaceLab", WorldLog, false, 0);
//         }
        //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    }

    // Return WorldPhys
    G4VPhysicalVolume* WorldPhys = new G4PVPlacement(NULL, G4ThreeVector(), WorldLog, "World", 0, false, 0);
    return WorldPhys;

}

void POLAR2DetectorConstruction::ConstructSDandField() {
//     POLAR2GlobalConfig* fPOLAR2GlobalConfig = POLAR2GlobalConfig::Instance();
//     if (!fPOLAR2GlobalConfig->primary_only) {
        POLAR2SensitiveDetector* POLAR2_SD = new POLAR2SensitiveDetector("POLAR2SD", "POLAR2HitsCollection");
        G4SDManager::GetSDMpointer()->AddNewDetector(POLAR2_SD);
        SetSensitiveDetector(SiPMLog_, POLAR2_SD);
        SetSensitiveDetector(LEAPPMTLog_, POLAR2_SD);
        //SetSensitiveDetector(ScintBlockLog_, POLAR2_SD);  // test the scintillation efficiency
//     }
}

void POLAR2DetectorConstruction::DefineMaterials_() {
    //-----------------------------------------------------------------------------------------------------------
    //--------Define Material
    //-----------------------------------------------------------------------------------------------------------

    POLAR2GlobalConfig* fPOLAR2GlobalConfig = POLAR2GlobalConfig::Instance();

    G4NistManager* man = G4NistManager::Instance();

    G4double density;                         // matter density
    G4double fractionmass;                    // compound proportion
    G4String name;                            // material name
    G4int ncomponents;                        // component number
    G4int nelements;
    G4int natoms;                             // atom proportion

    //Environment_
//     Environment_ = man->FindOrBuildMaterial("G4_Galactic");  // for space (Vacuum)
    Environment_ = man->FindOrBuildMaterial("G4_AIR");       // for ground (Air)

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    G4Element* elB  = man->FindOrBuildElement("B");
    G4Element* elC  = man->FindOrBuildElement("C");
    G4Element* elH  = man->FindOrBuildElement("H");
    G4Element* elHe  = man->FindOrBuildElement("He");
    G4Element* elO  = man->FindOrBuildElement("O");
    G4Element* elF  = man->FindOrBuildElement("F");
    G4Element* elK  = man->FindOrBuildElement("K");
    G4Element* elN  = man->FindOrBuildElement("N");
    G4Element* elAl = man->FindOrBuildElement("Al");
    G4Element* elAr = man->FindOrBuildElement("Ar");
    G4Element* elNa = man->FindOrBuildElement("Na");
    G4Element* elSi = man->FindOrBuildElement("Si");
    G4Element* elFe = man->FindOrBuildElement("Fe");
    G4Element* elCr = man->FindOrBuildElement("Cr");
    G4Element* elMn = man->FindOrBuildElement("Mn");
    G4Element* elMg = man->FindOrBuildElement("Mg");
    G4Element* elCs = man->FindOrBuildElement("Cs");
    G4Element* elSb = man->FindOrBuildElement("Sb");
    G4Element* elCu = man->FindOrBuildElement("Cu");
    G4Element* elZn = man->FindOrBuildElement("Zn");
    G4Element* elTi = man->FindOrBuildElement("Ti");
    G4Element* elZr = man->FindOrBuildElement("Zr");
    G4Element* elV  = man->FindOrBuildElement("V");

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    // theScint
    theScint = new G4Material(name = "EJ_248", density = 1.023*g/cm3, ncomponents = 2);
    theScint->AddElement(elH, fractionmass = 0.08482);  // 5.18E22 atoms per cm3
    theScint->AddElement(elC, fractionmass = 0.91518);  // 4.69E22 atoms per cm3
    theScint->GetIonisation()->SetBirksConstant(fPOLAR2GlobalConfig->birks_constant*mm/MeV);
    if (use_EJ200){
      theScint = new G4Material(name = "EJ_200", density = 1.023*g/cm3, ncomponents = 2);
      theScint->AddElement(elH, fractionmass = 0.08467);  // 5.17E22 atoms per cm3
      theScint->AddElement(elC, fractionmass = 0.91533);  // 4.69E22 atoms per cm3
      theScint->GetIonisation()->SetBirksConstant(fPOLAR2GlobalConfig->birks_constant*mm/MeV);
    }
        //FR4 for PCB
    //Silicon Oxide
    G4Material* SiliconOxide = new G4Material("SiliconOxide", density = 2.65*g/cm3, ncomponents = 2);
    SiliconOxide->AddElement(elSi, natoms=1);
    SiliconOxide->AddElement(elO,  natoms=2);
    //Diglycidyl Ether of Bisphenol A (First compound of epoxy resin Epotek 301-1)
    G4Material* Epoxy_1 = new G4Material("Epoxy_1", density = 1.16*g/cm3, ncomponents = 3);
    Epoxy_1->AddElement(elC, natoms=19);
    Epoxy_1->AddElement(elH, natoms=20);
    Epoxy_1->AddElement(elO, natoms=4);
    //1,4-Butanediol Diglycidyl Ether (Second compound of epoxy resin Epotek 301-1)
    G4Material* Epoxy_2 = new G4Material("Epoxy_2", density = 1.10*g/cm3, ncomponents = 3);
    Epoxy_2->AddElement(elC, natoms=10);
    Epoxy_2->AddElement(elH, natoms=18);
    Epoxy_2->AddElement(elO, natoms=4);
    //1,6-Hexanediamine 2,2,4-trimetyl (Third compound of epoxy resin Epotek 301-1)
    G4Material* Epoxy_3 = new G4Material("Epoxy_3", density = 1.16*g/cm3, ncomponents = 3);
    Epoxy_3->AddElement(elC, natoms=9);
    Epoxy_3->AddElement(elH, natoms=22);
    Epoxy_3->AddElement(elN, natoms=2);
    //Epoxy resin Epotek 301-1
    G4Material* Epoxy_Resin = new G4Material("Epoxy_Resin", density = 1.19*g/cm3, ncomponents = 3);
    Epoxy_Resin->AddMaterial(Epoxy_1, fractionmass=56*perCent);
    Epoxy_Resin->AddMaterial(Epoxy_2, fractionmass=24*perCent);
    Epoxy_Resin->AddMaterial(Epoxy_3, fractionmass=20*perCent);
    //FR4 PCB material
    FR4_ = new G4Material("FR4", density = 1.8*g/cm3, ncomponents=2);
    FR4_->AddMaterial(SiliconOxide, fractionmass=60*perCent);
    FR4_->AddMaterial(Epoxy_Resin,  fractionmass=40*perCent);

        // Accura25 (for the scintillator alignment grids), see https://pubchem.ncbi.nlm.nih.gov/ for compositions of molecules
    // 3,4- Epoxycyclohexylmethyl 3’,4’- epoxycyclohexane carboxylate 
    G4Material* Accura_1 = new G4Material("Accura_1", density = 1.17*g/cm3, ncomponents = 3);
    Accura_1->AddElement(elC, natoms=14);
    Accura_1->AddElement(elH, natoms=20);
    Accura_1->AddElement(elO, natoms=4);
    // 1,6-bis(2,3-epoxypropoxy)hexane 
    G4Material* Accura_2 = new G4Material("Accura_2", density = 1.076*g/cm3, ncomponents = 3);
    Accura_2->AddElement(elC, natoms=12);
    Accura_2->AddElement(elH, natoms=22);
    Accura_2->AddElement(elO, natoms=4);
    // Accura25 grid material
    Accura25 = new G4Material("Accura25", 1.19*g/cm3, ncomponents=2);
    Accura25->AddMaterial(Accura_1, fractionmass=60*perCent); //30-50%
    Accura25->AddMaterial(Accura_2, fractionmass=40*perCent); //15-30% + 2-4% of triarylsulfonium salts mixed with propylene carbonate, the rest is unknown

    // DC93-500 Dow Corning (for optical coupling)
    dc93_500 = new G4Material("Optical_Pad", 1.08*g/cm3, nelements=2);
    dc93_500->AddElement(elSi, natoms=1);
    dc93_500->AddElement(elO, natoms=2);
  
    // Epoxy Resin of the SiPM entrance window (S13361-6075NE-04), no composition provided, only density and refrindex, using the same epoxy resin used in FR4
    EpoxyResin = new G4Material("EpoxyResin", 1.21*g/cm3, ncomponents=3);
    EpoxyResin->AddMaterial(Epoxy_1, fractionmass=56*perCent);
    EpoxyResin->AddMaterial(Epoxy_2, fractionmass=24*perCent);
    EpoxyResin->AddMaterial(Epoxy_3, fractionmass=20*perCent);
  
    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......   

    const G4int n3Entries = 3;
    G4double Photon3Energy[n3Entries] = {2.48*eV, 2.93*eV, 3.10*eV};
    const G4int n25Entries = 25;
    G4double Photon25Energy[n25Entries] = {2.50*eV, 2.53*eV, 2.57*eV, 2.62*eV, 2.67*eV, 
									 2.72*eV, 2.75*eV, 2.77*eV, 2.79*eV, 2.81*eV,
									 2.84*eV, 2.87*eV, 2.88*eV, 2.90*eV, 2.91*eV, 
									 2.93*eV, 2.94*eV, 2.96*eV, 2.98*eV, 3.00*eV,
									 3.00*eV, 3.02*eV, 3.04*eV, 3.07*eV, 3.09*eV};
    
    //Set Optical Properties
      // Scintillator bar
    G4double PhotonEnergyEJ200[200] = {2.486*eV, 2.490*eV, 2.494*eV, 2.498*eV, 2.502*eV, 2.505*eV, 2.509*eV, 2.513*eV, 2.517*eV, 2.521*eV, 2.525*eV, 2.529*eV, 2.533*eV, 2.537*eV, 2.540*eV, 2.544*eV, 2.548*eV, 2.552*eV, 2.556*eV, 2.560*eV, 2.564*eV, 2.568*eV, 2.572*eV, 2.576*eV, 2.579*eV, 2.583*eV, 2.587*eV, 2.591*eV, 2.595*eV, 2.599*eV, 2.603*eV, 2.607*eV, 2.611*eV, 2.614*eV, 2.618*eV, 2.622*eV, 2.626*eV, 2.630*eV, 2.634*eV, 2.638*eV, 2.642*eV, 2.646*eV, 2.649*eV, 2.653*eV, 2.657*eV, 2.661*eV, 2.665*eV, 2.669*eV, 2.673*eV, 2.677*eV, 2.681*eV, 2.685*eV, 2.688*eV, 2.692*eV, 2.696*eV, 2.700*eV, 2.704*eV, 2.708*eV, 2.712*eV, 2.716*eV, 2.720*eV, 2.723*eV, 2.727*eV, 2.731*eV, 2.735*eV, 2.739*eV, 2.743*eV, 2.747*eV, 2.751*eV, 2.755*eV, 2.758*eV, 2.762*eV, 2.766*eV, 2.770*eV, 2.774*eV, 2.778*eV, 2.782*eV, 2.786*eV, 2.790*eV, 2.794*eV, 2.797*eV, 2.801*eV, 2.805*eV, 2.809*eV, 2.813*eV, 2.817*eV, 2.821*eV, 2.825*eV, 2.829*eV, 2.832*eV, 2.836*eV, 2.840*eV, 2.844*eV, 2.848*eV, 2.852*eV, 2.856*eV, 2.860*eV, 2.864*eV, 2.867*eV, 2.871*eV, 2.875*eV, 2.879*eV, 2.883*eV, 2.887*eV, 2.891*eV, 2.895*eV, 2.899*eV, 2.903*eV, 2.906*eV, 2.910*eV, 2.914*eV, 2.918*eV, 2.922*eV, 2.926*eV, 2.930*eV, 2.934*eV, 2.938*eV, 2.941*eV, 2.945*eV, 2.949*eV, 2.953*eV, 2.957*eV, 2.961*eV, 2.965*eV, 2.969*eV, 2.973*eV, 2.977*eV, 2.980*eV, 2.984*eV, 2.988*eV, 2.992*eV, 2.996*eV, 3.000*eV, 3.004*eV, 3.008*eV, 3.012*eV, 3.015*eV, 3.019*eV, 3.023*eV, 3.027*eV, 3.031*eV, 3.035*eV, 3.039*eV, 3.043*eV, 3.047*eV, 3.050*eV, 3.054*eV, 3.058*eV, 3.062*eV, 3.066*eV, 3.070*eV, 3.074*eV, 3.078*eV, 3.082*eV, 3.086*eV, 3.089*eV, 3.093*eV, 3.097*eV, 3.101*eV, 3.105*eV, 3.109*eV, 3.113*eV, 3.117*eV, 3.121*eV, 3.124*eV, 3.128*eV, 3.132*eV, 3.136*eV, 3.140*eV, 3.144*eV, 3.148*eV, 3.152*eV, 3.156*eV, 3.159*eV, 3.163*eV, 3.167*eV, 3.171*eV, 3.175*eV, 3.179*eV, 3.183*eV, 3.187*eV, 3.191*eV, 3.195*eV, 3.198*eV, 3.202*eV, 3.206*eV, 3.210*eV, 3.214*eV, 3.218*eV, 3.222*eV, 3.226*eV, 3.230*eV, 3.233*eV, 3.237*eV, 3.241*eV, 3.245*eV, 3.249*eV, 3.253*eV, 3.257*eV, 3.261*eV};
    G4double PhotonEnergyEJ248M[200] = {2.388*eV, 2.392*eV, 2.395*eV, 2.399*eV, 2.402*eV, 2.406*eV, 2.409*eV, 2.413*eV, 2.417*eV, 2.420*eV, 2.424*eV, 2.427*eV, 2.431*eV, 2.434*eV, 2.438*eV, 2.442*eV, 2.445*eV, 2.449*eV, 2.452*eV, 2.456*eV, 2.459*eV, 2.463*eV, 2.466*eV, 2.470*eV, 2.474*eV, 2.477*eV, 2.481*eV, 2.484*eV, 2.488*eV, 2.491*eV, 2.495*eV, 2.499*eV, 2.502*eV, 2.506*eV, 2.509*eV, 2.513*eV, 2.516*eV, 2.520*eV, 2.523*eV, 2.527*eV, 2.531*eV, 2.534*eV, 2.538*eV, 2.541*eV, 2.545*eV, 2.548*eV, 2.552*eV, 2.555*eV, 2.559*eV, 2.563*eV, 2.566*eV, 2.570*eV, 2.573*eV, 2.577*eV, 2.580*eV, 2.584*eV, 2.588*eV, 2.591*eV, 2.595*eV, 2.598*eV, 2.602*eV, 2.605*eV, 2.609*eV, 2.612*eV, 2.616*eV, 2.620*eV, 2.623*eV, 2.627*eV, 2.630*eV, 2.634*eV, 2.637*eV, 2.641*eV, 2.645*eV, 2.648*eV, 2.652*eV, 2.655*eV, 2.659*eV, 2.662*eV, 2.666*eV, 2.669*eV, 2.673*eV, 2.677*eV, 2.680*eV, 2.684*eV, 2.687*eV, 2.691*eV, 2.694*eV, 2.698*eV, 2.702*eV, 2.705*eV, 2.709*eV, 2.712*eV, 2.716*eV, 2.719*eV, 2.723*eV, 2.726*eV, 2.730*eV, 2.734*eV, 2.737*eV, 2.741*eV, 2.744*eV, 2.748*eV, 2.751*eV, 2.755*eV, 2.758*eV, 2.762*eV, 2.766*eV, 2.769*eV, 2.773*eV, 2.776*eV, 2.780*eV, 2.783*eV, 2.787*eV, 2.791*eV, 2.794*eV, 2.798*eV, 2.801*eV, 2.805*eV, 2.808*eV, 2.812*eV, 2.815*eV, 2.819*eV, 2.823*eV, 2.826*eV, 2.830*eV, 2.833*eV, 2.837*eV, 2.840*eV, 2.844*eV, 2.848*eV, 2.851*eV, 2.855*eV, 2.858*eV, 2.862*eV, 2.865*eV, 2.869*eV, 2.872*eV, 2.876*eV, 2.880*eV, 2.883*eV, 2.887*eV, 2.890*eV, 2.894*eV, 2.897*eV, 2.901*eV, 2.905*eV, 2.908*eV, 2.912*eV, 2.915*eV, 2.919*eV, 2.922*eV, 2.926*eV, 2.929*eV, 2.933*eV, 2.937*eV, 2.940*eV, 2.944*eV, 2.947*eV, 2.951*eV, 2.954*eV, 2.958*eV, 2.962*eV, 2.965*eV, 2.969*eV, 2.972*eV, 2.976*eV, 2.979*eV, 2.983*eV, 2.986*eV, 2.990*eV, 2.994*eV, 2.997*eV, 3.001*eV, 3.004*eV, 3.008*eV, 3.011*eV, 3.015*eV, 3.018*eV, 3.022*eV, 3.026*eV, 3.029*eV, 3.033*eV, 3.036*eV, 3.040*eV, 3.043*eV, 3.047*eV, 3.051*eV, 3.054*eV, 3.058*eV, 3.061*eV, 3.065*eV, 3.068*eV, 3.072*eV, 3.075*eV, 3.079*eV, 3.083*eV, 3.086*eV, 3.090*eV, 3.093*eV, 3.097*eV};
    G4double Scintil_EJ248M[200] = {0.005, 0.005, 0.004, 0.005, 0.010, 0.013, 0.008, 0.011, 0.016, 0.015, 0.022, 0.024, 0.026, 0.025, 0.027, 0.032, 0.029, 0.030, 0.033, 0.034, 0.036, 0.041, 0.045, 0.050, 0.053, 0.056, 0.059, 0.063, 0.066, 0.068, 0.075, 0.075, 0.075, 0.078, 0.084, 0.088, 0.092, 0.096, 0.100, 0.100, 0.105, 0.110, 0.115, 0.118, 0.123, 0.130, 0.134, 0.139, 0.144, 0.150, 0.152, 0.157, 0.165, 0.171, 0.176, 0.180, 0.190, 0.202, 0.208, 0.216, 0.221, 0.227, 0.237, 0.247, 0.257, 0.261, 0.270, 0.283, 0.294, 0.304, 0.313, 0.328, 0.335, 0.348, 0.356, 0.367, 0.375, 0.383, 0.390, 0.401, 0.410, 0.416, 0.420, 0.424, 0.425, 0.430, 0.434, 0.438, 0.444, 0.453, 0.455, 0.457, 0.466, 0.468, 0.468, 0.477, 0.483, 0.491, 0.502, 0.509, 0.512, 0.517, 0.523, 0.528, 0.536, 0.546, 0.555, 0.562, 0.573, 0.582, 0.591, 0.601, 0.608, 0.614, 0.634, 0.640, 0.649, 0.658, 0.672, 0.683, 0.693, 0.704, 0.716, 0.728, 0.736, 0.750, 0.758, 0.769, 0.786, 0.792, 0.796, 0.803, 0.810, 0.819, 0.824, 0.829, 0.832, 0.843, 0.854, 0.855, 0.857, 0.862, 0.867, 0.873, 0.878, 0.880, 0.879, 0.877, 0.875, 0.873, 0.866, 0.861, 0.855, 0.843, 0.833, 0.815, 0.801, 0.781, 0.750, 0.722, 0.696, 0.661, 0.630, 0.607, 0.577, 0.530, 0.514, 0.478, 0.457, 0.424, 0.402, 0.373, 0.341, 0.313, 0.297, 0.273, 0.255, 0.238, 0.221, 0.199, 0.182, 0.169, 0.151, 0.138, 0.121, 0.105, 0.095, 0.085, 0.069, 0.058, 0.050, 0.040, 0.031, 0.024, 0.019, 0.014, 0.012, 0.007, 0.003, 0.005};
    G4double Scintil_EJ200[200] = {0.071, 0.071, 0.072, 0.074, 0.080, 0.086, 0.089, 0.091, 0.096, 0.103, 0.107, 0.110, 0.113, 0.118, 0.123, 0.126, 0.132, 0.138, 0.144, 0.152, 0.158, 0.163, 0.171, 0.177, 0.184, 0.189, 0.197, 0.205, 0.212, 0.219, 0.227, 0.235, 0.249, 0.258, 0.266, 0.273, 0.286, 0.298, 0.311, 0.322, 0.335, 0.347, 0.364, 0.375, 0.382, 0.397, 0.405, 0.416, 0.426, 0.435, 0.441, 0.447, 0.454, 0.461, 0.464, 0.466, 0.472, 0.476, 0.478, 0.487, 0.498, 0.503, 0.506, 0.512, 0.519, 0.529, 0.540, 0.548, 0.554, 0.561, 0.569, 0.577, 0.588, 0.601, 0.608, 0.616, 0.624, 0.636, 0.651, 0.659, 0.674, 0.682, 0.699, 0.713, 0.728, 0.737, 0.748, 0.760, 0.772, 0.787, 0.801, 0.812, 0.824, 0.837, 0.848, 0.858, 0.868, 0.878, 0.885, 0.892, 0.898, 0.905, 0.911, 0.918, 0.922, 0.929, 0.935, 0.940, 0.945, 0.945, 0.944, 0.942, 0.940, 0.935, 0.928, 0.918, 0.908, 0.891, 0.871, 0.847, 0.814, 0.788, 0.748, 0.715, 0.676, 0.640, 0.605, 0.567, 0.527, 0.491, 0.465, 0.426, 0.396, 0.366, 0.337, 0.311, 0.283, 0.264, 0.245, 0.220, 0.202, 0.185, 0.172, 0.149, 0.131, 0.113, 0.102, 0.083, 0.067, 0.056, 0.040, 0.032, 0.025, 0.018, 0.012, 0.010, 0.009, 0.007, 0.005, 0.004, 0.004, 0.003, 0.003, 0.002, 0.002, 0.002, 0.002, 0.002, 0.000, -0.001, -0.000, 0.002, 0.001, 0.000, 0.001, 0.002, 0.002, 0.001, 0.001, 0.000, 0.001, 0.001, 0.000, 0.000, 0.000, 0.001, 0.003, 0.003, 0.003, 0.002, 0.001, 0.000, 0.001, 0.002, 0.002, 0.001, 0.000, 0.001, 0.002, 0.002};

    G4double SciRefrInd[n25Entries] =	{1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58,
									 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58,
									 1.58, 1.58, 1.58, 1.58, 1.58};

    G4double SciAbs_EJ248M[n25Entries] =	{250.0*cm, 250.0*cm, 250.0*cm, 250.0*cm, 250.0*cm, 
                    250.0*cm, 250.0*cm, 250.0*cm, 250.0*cm, 250.0*cm, 
                    250.0*cm, 250.0*cm, 250.0*cm, 250.0*cm, 250.0*cm,
                    250.0*cm, 250.0*cm, 250.0*cm, 250.0*cm, 250.0*cm,
                    250.0*cm, 250.0*cm, 250.0*cm, 250.0*cm, 250.0*cm};

    // test the scintillation efficiency -> Scint has sensitive volume
    G4double SciAbs_sensitive[n25Entries] =	{0.001*mm, 0.001*mm, 0.001*mm, 0.001*mm, 0.001*mm,
                    0.001*mm, 0.001*mm, 0.001*mm, 0.001*mm, 0.001*mm,
                    0.001*mm, 0.001*mm, 0.001*mm, 0.001*mm, 0.001*mm,
                    0.001*mm, 0.001*mm, 0.001*mm, 0.001*mm, 0.001*mm,
                    0.001*mm, 0.001*mm, 0.001*mm, 0.001*mm, 0.001*mm};

    G4double SciAbs_EJ200[n25Entries] =	{380.0*cm, 380.0*cm, 380.0*cm, 380.0*cm, 380.0*cm, 
									 380.0*cm, 380.0*cm, 380.0*cm, 380.0*cm, 380.0*cm, 
									 380.0*cm, 380.0*cm, 380.0*cm, 380.0*cm, 380.0*cm,
									 380.0*cm, 380.0*cm, 380.0*cm, 380.0*cm, 380.0*cm,
									 380.0*cm, 380.0*cm, 380.0*cm, 380.0*cm, 380.0*cm};


    SciMPT = new G4MaterialPropertiesTable();
    SciMPT->AddProperty("RINDEX",       Photon25Energy, SciRefrInd,  n25Entries);

    //SciMPT->AddProperty("ABSLENGTH",    Photon25Energy, SciAbs_sensitive, n25Entries);
    
    if (use_EJ200){
      SciMPT->AddProperty("ABSLENGTH",    Photon25Energy, SciAbs_EJ200, n25Entries);
      SciMPT->AddProperty("FASTCOMPONENT",PhotonEnergyEJ200, Scintil_EJ200, 200);
      SciMPT->AddConstProperty("SCINTILLATIONYIELD",10000./MeV);
    }else{
      SciMPT->AddProperty("ABSLENGTH",    Photon25Energy, SciAbs_EJ248M, n25Entries);
      SciMPT->AddProperty("FASTCOMPONENT",PhotonEnergyEJ248M, Scintil_EJ248M, 200);
      SciMPT->AddConstProperty("SCINTILLATIONYIELD",9200./MeV);		//1./(90.0*eV));	//,1./keV);
    }

    SciMPT->AddConstProperty("RESOLUTIONSCALE",1.0);
    SciMPT->AddConstProperty("FASTSCINTILLATIONRISETIME", 0.9*ns); //rise time of the scintillators
    SciMPT->AddConstProperty("FASTTIMECONSTANT", 2.1*ns); //decay time of the scintillators (only fast component in our case, no slow component)
    SciMPT->AddConstProperty("YIELDRATIO",1.0); // fast/(fast+slow)
  
    theScint -> SetMaterialPropertiesTable(SciMPT);

 
      // Air MPT
    G4double AirRefrInd[n3Entries] = { 1.00, 1.00, 1.00};
    G4MaterialPropertiesTable* AirMPT = new G4MaterialPropertiesTable();
    AirMPT->AddProperty("RINDEX", Photon3Energy, AirRefrInd, n3Entries);
    Environment_->SetMaterialPropertiesTable(AirMPT); 
    
      //Glass MPT
    G4double SiPM_RIND[n3Entries]={1.51, 1.51, 1.51};	//Ref Index between 1.51 - 1.54 according to wikipedia  
    G4double SiPMAbs[n3Entries] =	{0.001*mm, 0.001*mm, 0.001*mm};
    //G4double SiPM_EFF[n3Entries] = {1.,1.,1.};
    G4MaterialPropertiesTable *SiPMMPT = new G4MaterialPropertiesTable();
    SiPMMPT->AddProperty("RINDEX", Photon3Energy, SiPM_RIND,n3Entries);
    SiPMMPT->AddProperty("ABSLENGTH", Photon3Energy, SiPMAbs,n3Entries);
    //SiPMMPT->AddProperty("EFFICIENCY", Photon3Energy, SiPM_EFF,n3Entries);

    FR4_->SetMaterialPropertiesTable(SiPMMPT);
    
      // SiPM Surface (for PDE)
    OpResinSiPMSurface = new G4OpticalSurface("ResinSiPMSurface");  
    OpResinSiPMSurface -> SetType(dielectric_metal); 
    OpResinSiPMSurface -> SetFinish(polished);
    OpResinSiPMSurface -> SetModel(unified);

    // SiPM: 0% reflectivity, efficiency=absorbed fraction=PDE from Hamamatsu datasheet (S13360-6075)

    G4double reflSiPM[n3Entries] = {0.0, 0.0, 0.0};   

    G4double effSiPM_photonenergy[200] = {1.374*eV, 1.386*eV, 1.399*eV, 1.411*eV, 1.424*eV, 1.436*eV, 1.449*eV, 1.461*eV, 1.474*eV, 1.486*eV, 1.499*eV, 1.511*eV, 1.524*eV, 1.536*eV, 1.549*eV, 1.562*eV, 1.574*eV, 1.587*eV, 1.599*eV, 1.612*eV, 1.624*eV, 1.637*eV, 1.649*eV, 1.662*eV, 1.674*eV, 1.687*eV, 1.699*eV, 1.712*eV, 1.724*eV, 1.737*eV, 1.749*eV, 1.762*eV, 1.774*eV, 1.787*eV, 1.799*eV, 1.812*eV, 1.824*eV, 1.837*eV, 1.849*eV, 1.862*eV, 1.874*eV, 1.887*eV, 1.900*eV, 1.912*eV, 1.925*eV, 1.937*eV, 1.950*eV, 1.962*eV, 1.975*eV, 1.987*eV, 2.000*eV, 2.012*eV, 2.025*eV, 2.037*eV, 2.050*eV, 2.062*eV, 2.075*eV, 2.087*eV, 2.100*eV, 2.112*eV, 2.125*eV, 2.137*eV, 2.150*eV, 2.162*eV, 2.175*eV, 2.187*eV, 2.200*eV, 2.213*eV, 2.225*eV, 2.238*eV, 2.250*eV, 2.263*eV, 2.275*eV, 2.288*eV, 2.300*eV, 2.313*eV, 2.325*eV, 2.338*eV, 2.350*eV, 2.363*eV, 2.375*eV, 2.388*eV, 2.400*eV, 2.413*eV, 2.425*eV, 2.438*eV, 2.450*eV, 2.463*eV, 2.475*eV, 2.488*eV, 2.500*eV, 2.513*eV, 2.525*eV, 2.538*eV, 2.551*eV, 2.563*eV, 2.576*eV, 2.588*eV, 2.601*eV, 2.613*eV, 2.626*eV, 2.638*eV, 2.651*eV, 2.663*eV, 2.676*eV, 2.688*eV, 2.701*eV, 2.713*eV, 2.726*eV, 2.738*eV, 2.751*eV, 2.763*eV, 2.776*eV, 2.788*eV, 2.801*eV, 2.813*eV, 2.826*eV, 2.838*eV, 2.851*eV, 2.863*eV, 2.876*eV, 2.889*eV, 2.901*eV, 2.914*eV, 2.926*eV, 2.939*eV, 2.951*eV, 2.964*eV, 2.976*eV, 2.989*eV, 3.001*eV, 3.014*eV, 3.026*eV, 3.039*eV, 3.051*eV, 3.064*eV, 3.076*eV, 3.089*eV, 3.101*eV, 3.114*eV, 3.126*eV, 3.139*eV, 3.151*eV, 3.164*eV, 3.176*eV, 3.189*eV, 3.201*eV, 3.214*eV, 3.227*eV, 3.239*eV, 3.252*eV, 3.264*eV, 3.277*eV, 3.289*eV, 3.302*eV, 3.314*eV, 3.327*eV, 3.339*eV, 3.352*eV, 3.364*eV, 3.377*eV, 3.389*eV, 3.402*eV, 3.414*eV, 3.427*eV, 3.439*eV, 3.452*eV, 3.464*eV, 3.477*eV, 3.489*eV, 3.502*eV, 3.514*eV, 3.527*eV, 3.540*eV, 3.552*eV, 3.565*eV, 3.577*eV, 3.590*eV, 3.602*eV, 3.615*eV, 3.627*eV, 3.640*eV, 3.652*eV, 3.665*eV, 3.677*eV, 3.690*eV, 3.702*eV, 3.715*eV, 3.727*eV, 3.740*eV, 3.752*eV, 3.765*eV, 3.777*eV, 3.790*eV, 3.802*eV, 3.815*eV, 3.827*eV, 3.840*eV, 3.852*eV, 3.865*eV};    //G4double effSiPM[178] = {26.809, 28.931, 31.150, 32.517, 35.069, 38.270, 39.929, 42.432, 44.545, 46.076, 48.345, 50.629, 53.858, 55.369, 56.041, 58.169, 60.373, 62.452, 63.548, 65.517, 67.175, 69.871, 71.144, 72.425, 72.731, 73.211, 74.893, 75.976, 77.003, 77.973, 79.189, 79.453, 80.888, 81.685, 82.448, 83.458, 83.817, 84.744, 85.018, 87.197, 85.597, 87.987, 88.895, 90.738, 92.032, 92.774, 93.351, 93.892, 93.997, 94.361, 94.773, 95.289, 95.722, 96.230, 97.176, 97.635, 98.056, 98.227, 98.366, 99.588, 99.259, 99.059, 99.744, 99.934, 100.367, 100.391, 100.108, 99.792, 99.939, 100.251, 99.381, 99.648, 99.543, 98.951, 98.314, 98.661, 97.997, 96.654, 95.963, 94.603, 93.270, 91.565, 88.877, 88.192, 86.680, 83.807, 82.255, 80.589, 78.944, 76.490, 74.652, 72.964, 70.634, 68.774, 66.947, 65.559, 63.678, 62.226, 60.026, 58.630, 57.835, 56.575, 55.231, 54.009, 52.955, 52.109, 50.834, 49.814, 48.709, 47.998, 46.531, 46.880, 45.744, 44.870, 44.098, 43.241, 42.392, 41.461, 41.095, 40.820, 39.898, 39.173, 38.814, 37.743, 36.603, 35.474, 34.373, 33.275, 32.129, 31.540, 30.425, 29.629, 28.760, 27.820, 26.855, 26.102, 25.546, 24.885, 24.491, 23.352, 22.964, 22.344, 21.801, 20.864, 20.253, 19.569, 19.038, 18.292, 17.840, 17.389, 16.862, 16.260, 15.389, 14.791, 13.684, 13.316, 12.630, 12.098, 11.446, 10.939, 10.287, 10.566, 9.414, 8.846, 8.563, 7.864, 7.454, 7.036, 6.409, 5.874, 5.211, 4.761, 4.386, 3.846, 3.171, 2.834, 3.414, 2.457};
    G4double effSiPM[200] = {0.115, 0.125, 0.136, 0.141, 0.148, 0.160, 0.169, 0.178, 0.185, 0.191, 0.196, 0.201, 0.208, 0.212, 0.216, 0.227, 0.235, 0.238, 0.242, 0.250, 0.257, 0.264, 0.269, 0.275, 0.283, 0.290, 0.298, 0.302, 0.305, 0.310, 0.311, 0.314, 0.323, 0.327, 0.331, 0.336, 0.340, 0.341, 0.348, 0.352, 0.357, 0.359, 0.362, 0.364, 0.366, 0.376, 0.379, 0.386, 0.395, 0.397, 0.399, 0.402, 0.402, 0.403, 0.404, 0.405, 0.408, 0.409, 0.411, 0.415, 0.417, 0.418, 0.420, 0.421, 0.422, 0.426, 0.425, 0.424, 0.427, 0.428, 0.430, 0.430, 0.429, 0.429, 0.428, 0.428, 0.429, 0.428, 0.426, 0.427, 0.426, 0.425, 0.424, 0.422, 0.422, 0.420, 0.418, 0.415, 0.413, 0.411, 0.408, 0.405, 0.401, 0.398, 0.395, 0.393, 0.390, 0.387, 0.384, 0.381, 0.379, 0.376, 0.373, 0.370, 0.367, 0.363, 0.360, 0.356, 0.353, 0.349, 0.345, 0.342, 0.338, 0.334, 0.330, 0.326, 0.322, 0.318, 0.315, 0.311, 0.307, 0.304, 0.300, 0.296, 0.293, 0.289, 0.285, 0.281, 0.277, 0.274, 0.270, 0.266, 0.263, 0.260, 0.257, 0.253, 0.249, 0.245, 0.242, 0.238, 0.235, 0.232, 0.227, 0.222, 0.219, 0.215, 0.210, 0.206, 0.199, 0.195, 0.191, 0.187, 0.183, 0.177, 0.175, 0.171, 0.166, 0.161, 0.159, 0.156, 0.151, 0.149, 0.145, 0.141, 0.137, 0.136, 0.132, 0.129, 0.125, 0.122, 0.118, 0.116, 0.114, 0.109, 0.104, 0.098, 0.095, 0.092, 0.089, 0.084, 0.076, 0.068, 0.064, 0.060, 0.056, 0.053, 0.044, 0.045, 0.041, 0.039, 0.035, 0.033, 0.027, 0.025, 0.023, 0.021, 0.018, 0.014, 0.014, 0.011};    //G4double SiPM_EFF[9] = {1.0,0.0,1.0,0.0,1.0,0.0,1.0,0.0,1.0};

    G4double perfectSiPM[n3Entries] = {1, 1, 1}; // 100% PDE, perfect SiPM

    G4MaterialPropertiesTable *ResinSiPMProperty = new G4MaterialPropertiesTable();
    ResinSiPMProperty -> AddProperty("REFLECTIVITY", Photon3Energy, reflSiPM, n3Entries);
    
    if (use_SiPM_PDE){
      ResinSiPMProperty -> AddProperty("EFFICIENCY", effSiPM_photonenergy, effSiPM, 200);
    }else{
      ResinSiPMProperty -> AddProperty("EFFICIENCY", Photon3Energy, perfectSiPM, n3Entries);
    }

    OpResinSiPMSurface -> SetMaterialPropertiesTable(ResinSiPMProperty);
     
      // SiPM to resin Surface (if photon going out of SiPM)
    OpSiPMResinSurface = new G4OpticalSurface("SiPMResinSurface");  
    OpSiPMResinSurface -> SetType(dielectric_metal); 
    OpSiPMResinSurface -> SetFinish(polished);
    OpSiPMResinSurface -> SetModel(unified);
    
    G4double reflSiPMResin[n3Entries] = {1.0, 1.0, 1.0};   

    G4MaterialPropertiesTable *SiPMResinProperty = new G4MaterialPropertiesTable();
    SiPMResinProperty -> AddProperty("REFLECTIVITY", Photon3Energy, reflSiPMResin, n3Entries);
    OpSiPMResinSurface -> SetMaterialPropertiesTable(SiPMResinProperty);

      // Optical pad to resin
    OpPadResinSurface = new G4OpticalSurface("PadResinSurface");  
    OpPadResinSurface -> SetType(dielectric_dielectric); 
    OpPadResinSurface -> SetFinish(ground);
    OpPadResinSurface -> SetModel(unified);

      // Resin to optical pad
    OpResinPadSurface = new G4OpticalSurface("ResinPadSurface");  
    OpResinPadSurface -> SetType(dielectric_dielectric); 
    OpResinPadSurface -> SetFinish(ground);
    OpResinPadSurface -> SetModel(unified);
    


    /*
    G4double density_foil = 2.7*g/cm3;    // Source; http://www.alufoil.org/eng/alufoil_345.html
    // Al Foil
    claryl = new G4Material("Aluminium_Foil", density_foil, 1);
    claryl->AddElement(elAl, natoms=1);
    
    reflcoating = new G4Material("Aluminium_Foil", density_foil, 1);
    reflcoating->AddElement(elAl, natoms=1);
    */


        // Vikuiti ESR 3M

    // Polyethylene naphthalate (poly(ethylene 2,6-naphthalate) or PEN)
    G4Material* ESR_PEN = new G4Material("PEN", density = 1.38*g/cm3, ncomponents = 3);
    ESR_PEN->AddElement(elC, natoms=14);
    ESR_PEN->AddElement(elH, natoms=10);
    ESR_PEN->AddElement(elO, natoms=4);
    // tetrafluoroethylene TFE
    G4Material* ESR_TFE = new G4Material("TFE", density = 1.519*g/cm3, ncomponents = 2);
    ESR_TFE->AddElement(elC, natoms=2);
    ESR_TFE->AddElement(elF, natoms=4);
    // hexafluoropropylene HFP
    G4Material* ESR_HFP = new G4Material("HFP", density = 1.33*g/cm3, ncomponents = 2);
    ESR_HFP->AddElement(elC, natoms=3);
    ESR_HFP->AddElement(elF, natoms=6);
    // vinylidene fluoride PVDF
    G4Material* ESR_PVDF = new G4Material("PVDF", density = 1.78*g/cm3, ncomponents = 3);
    ESR_PVDF->AddElement(elC, natoms=2);
    ESR_PVDF->AddElement(elF, natoms=2);  
    ESR_PVDF->AddElement(elH, natoms=2);  
    // 3M's Dyneon™ THV Fluorothermo-plastic
    G4Material* ESR_THV = new G4Material("THV", density = 1.97*g/cm3, ncomponents = 3);
    ESR_THV->AddMaterial(ESR_TFE, fractionmass=34*perCent);
    ESR_THV->AddMaterial(ESR_HFP, fractionmass=33*perCent);
    ESR_THV->AddMaterial(ESR_PVDF, fractionmass=33*perCent);
    // Vikuiti ESR (guessed as PEN and THV, see https://patents.google.com/patent/US20110228511), read also http://dx.doi.org/10.1016/j.nima.2017.01.051
    VikuitiESR_ = new G4Material("Vikuiti", 1.29*g/cm3, ncomponents=2);  // density confirmed with measurements
    VikuitiESR_->AddMaterial(ESR_PEN, fractionmass=50*perCent);
    VikuitiESR_->AddMaterial(ESR_THV, fractionmass=50*perCent);
    


    // PEEK_Plastic -> for LEAP housing
    PEEK_Plastic_ = new G4Material(name = "PEEK_Plastic", density = 1.32*g/cm3, ncomponents = 3);
    PEEK_Plastic_->AddElement(elC, natoms = 20);
    PEEK_Plastic_->AddElement(elH, natoms = 12);
    PEEK_Plastic_->AddElement(elO, natoms = 3);


    
      //Aluminum MPT
    G4double AluRefrInd[n3Entries] = { 1.39, 1.39, 1.39};  //Aluminum   
    // From: www.robinwood.com/Catalog/Technical/Gen3DTuts/Gen3DPages/RefractionIndexList.html
    G4MaterialPropertiesTable* AluMPT = new G4MaterialPropertiesTable();
    AluMPT->AddProperty("RINDEX", Photon3Energy, AluRefrInd, n3Entries);  
    //claryl->SetMaterialPropertiesTable(AluMPT);
    //reflcoating->SetMaterialPropertiesTable(AluMPT);
    VikuitiESR_->SetMaterialPropertiesTable(AluMPT);
    

    
      //Accura25 MPT
    G4double PP_RIND[n3Entries]={1.49, 1.49, 1.49};	//approximate value searched on the web
    G4double PP_Abs[n3Entries] = {1.166*mm, 1.166*mm, 1.166*mm};  // from transmittance measured at CERN

    G4MaterialPropertiesTable *PPMPT = new G4MaterialPropertiesTable();
    PPMPT->AddProperty("RINDEX", Photon3Energy, PP_RIND,n3Entries);
    PPMPT->AddProperty("ABSLENGTH", Photon3Energy, PP_Abs,n3Entries);
    Accura25->SetMaterialPropertiesTable(PPMPT);
    
    G4double dc93_500_RIND[n3Entries]={1.5324, 1.5324, 1.5324};	
    G4MaterialPropertiesTable *dc93_500MPT = new G4MaterialPropertiesTable();
    dc93_500MPT->AddProperty("RINDEX", Photon3Energy, dc93_500_RIND,n3Entries);
    dc93_500->SetMaterialPropertiesTable(dc93_500MPT);
    
      //Epoxy Resin MPT
    G4double EpoxyResin_RIND[n3Entries]={1.55, 1.55, 1.55};	//from Hamamatsu - 1.55 for epoxy s13360 (1.57 for s14160 and 1.41 for silicon s13360)
    G4MaterialPropertiesTable *EpoxyResinMPT = new G4MaterialPropertiesTable();
    EpoxyResinMPT->AddProperty("RINDEX", Photon3Energy, EpoxyResin_RIND,n3Entries);
    EpoxyResin->SetMaterialPropertiesTable(EpoxyResinMPT);
    
    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    // define other materials that are not sensitive


    // Air-Grid Surface
    OpAirGridSurface = new G4OpticalSurface("Air-GridSurface");  
    OpAirGridSurface -> SetType(dielectric_dielectric); 
    OpAirGridSurface -> SetFinish(polished);
    OpAirGridSurface -> SetModel(glisur);

    // Grid 60% transmittance at 0.2mm, reflectance probably low

    G4double notabsGrid[n3Entries] = {0.5, 0.5, 0.5};   //not absorbed fraction, fine tuning to get 60%
    G4double transGrid[n3Entries] = {1.0, 1.0, 1.0};   //transmitted fraction

    G4MaterialPropertiesTable *AirGridProperty = new G4MaterialPropertiesTable();
    AirGridProperty -> AddProperty("REFLECTIVITY", Photon3Energy, notabsGrid, n3Entries);
    AirGridProperty -> AddProperty("TRANSMITTANCE", Photon3Energy, transGrid, n3Entries);

    OpAirGridSurface -> SetMaterialPropertiesTable(AirGridProperty);
    
  //-------- Grid-Air Surface
    OpGridAirSurface = new G4OpticalSurface("Grid-AirSurface");  
    OpGridAirSurface -> SetType(dielectric_dielectric);
    OpGridAirSurface -> SetFinish(polished);
    OpGridAirSurface -> SetModel(glisur);
    //OpGridAirSurface -> SetFinish(ground);
    //G4double polish = 1.0;
    //OpGridAirSurface -> SetPolish(polish);

    G4double notabsGridAir[n3Entries] = {0.76, 0.76, 0.76}; // fine tuning to get 60%
    G4double transGridAir[n3Entries] = {1.0, 1.0, 1.0};

    G4MaterialPropertiesTable *GridAirProperty = new G4MaterialPropertiesTable();
    GridAirProperty -> AddProperty("REFLECTIVITY", Photon3Energy, notabsGridAir, n3Entries);
    //GridAirProperty -> AddProperty("TRANSMITTANCE", Photon3Energy, transGridAir, n3Entries);

    OpGridAirSurface -> SetMaterialPropertiesTable(GridAirProperty);


  // ############################ Perfect reflector for LY testing ############################
  // Air-PerfRefl Surface
    OpAirPerfReflSurface = new G4OpticalSurface("Air-PerfReflSurface");  
    OpAirPerfReflSurface -> SetType(dielectric_metal); 
    OpAirPerfReflSurface -> SetFinish(polished);
    OpAirPerfReflSurface -> SetModel(glisur);

    // Perfect reflector: 100% reflectivity, 0% transmittance

    G4double reflPerfRefl[n3Entries] = {1.0, 1.0, 1.0};   //not absorbed fraction = reflected (dielectric_metal)

    G4MaterialPropertiesTable *AirPerfReflProperty = new G4MaterialPropertiesTable();
    AirPerfReflProperty -> AddProperty("REFLECTIVITY", Photon3Energy, reflPerfRefl, n3Entries);

    OpAirPerfReflSurface -> SetMaterialPropertiesTable(AirPerfReflProperty);
     
 
  //-------- PerfRefl-Air Surface
    OpPerfReflAirSurface = new G4OpticalSurface("PerfRefl-AirSurface");  
    OpPerfReflAirSurface -> SetType(dielectric_metal);
    OpPerfReflAirSurface -> SetFinish(polished);
    OpPerfReflAirSurface -> SetModel(glisur);

    G4double reflAirPerfRefl[n3Entries] = {0.0, 0.0, 0.0};

    G4MaterialPropertiesTable *PerfReflAirProperty = new G4MaterialPropertiesTable();
    PerfReflAirProperty -> AddProperty("REFLECTIVITY", Photon3Energy, reflAirPerfRefl, n3Entries);

    OpPerfReflAirSurface -> SetMaterialPropertiesTable(PerfReflAirProperty);

  // ##########################################################################################


if (vikuiti_is_metal) {
  // ################ Vikuiti as dielectric_metal (with 1.5% absorbance, used for lightoutput analysis) ################
  // Air-Vikuiti Surface
    OpAirWrapSurface = new G4OpticalSurface("Air-VikuitiSurface");  
    OpAirWrapSurface -> SetType(dielectric_metal); 
    OpAirWrapSurface -> SetFinish(polished);
    OpAirWrapSurface -> SetModel(glisur);

    // Vikuiti: 98.5% reflectivity, 1.5% transmittance->here not transmittance but absorbance

    G4double reflVikuiti[n3Entries] = {0.985, 0.985, 0.985};   //not absorbed fraction = reflected (dielectric_metal)

    G4MaterialPropertiesTable *AirWrapProperty = new G4MaterialPropertiesTable();
    AirWrapProperty -> AddProperty("REFLECTIVITY", Photon3Energy, reflVikuiti, n3Entries);

    OpAirWrapSurface -> SetMaterialPropertiesTable(AirWrapProperty);
     
 
  //-------- PerfRefl-Air Surface
    OpWrapAirSurface = new G4OpticalSurface("Vikuiti-AirSurface");  
    OpWrapAirSurface -> SetType(dielectric_metal);
    OpWrapAirSurface -> SetFinish(polished);
    OpWrapAirSurface -> SetModel(glisur);

    G4double reflAirVikuiti[n3Entries] = {0.0, 0.0, 0.0};

    G4MaterialPropertiesTable *WrapAirProperty = new G4MaterialPropertiesTable();
    WrapAirProperty -> AddProperty("REFLECTIVITY", Photon3Energy, reflAirVikuiti, n3Entries);

    OpWrapAirSurface -> SetMaterialPropertiesTable(WrapAirProperty);

  // ##########################################################################################
}else{
  // ################ Vikuiti as dielectric_dielectric (with 1.5% transmittance, used for crosstalk analysis) ################
  // Air-Vikuiti Surface
    OpAirWrapSurface = new G4OpticalSurface("Air-VikuitiSurface");  
    OpAirWrapSurface -> SetType(dielectric_dielectric); 
    OpAirWrapSurface -> SetFinish(polished);
    OpAirWrapSurface -> SetModel(glisur);

    // Vikuiti: 98.5% reflectivity, 1.5% transmittance

    G4double notabsVikuiti[n3Entries] = {1.0, 1.0, 1.0};   //not absorbed fraction
    G4double transVikuiti[n3Entries] = {0.015, 0.015, 0.015};   //transmitted fraction

    G4MaterialPropertiesTable *AirWrapProperty = new G4MaterialPropertiesTable();
    AirWrapProperty -> AddProperty("REFLECTIVITY", Photon3Energy, notabsVikuiti, n3Entries);
    AirWrapProperty -> AddProperty("TRANSMITTANCE", Photon3Energy, transVikuiti, n3Entries);

    OpAirWrapSurface -> SetMaterialPropertiesTable(AirWrapProperty);
     
 
  //-------- Vikuiti-Air Surface
    OpWrapAirSurface = new G4OpticalSurface("Vikuiti-AirSurface");  
    OpWrapAirSurface -> SetType(dielectric_dielectric);
    OpWrapAirSurface -> SetFinish(polished);
    OpWrapAirSurface -> SetModel(glisur);
    //OpWrapAirSurface -> SetFinish(ground);
    //G4double polish = 1.0;
    //OpWrapAirSurface -> SetPolish(polish);

    G4double notabsAir[n3Entries] = {1.0, 1.0, 1.0};
    G4double transAir[n3Entries] = {1.0, 1.0, 1.0};

    G4MaterialPropertiesTable *WrapAirProperty = new G4MaterialPropertiesTable();
    WrapAirProperty -> AddProperty("REFLECTIVITY", Photon3Energy, notabsAir, n3Entries);
    WrapAirProperty -> AddProperty("TRANSMITTANCE", Photon3Energy, transAir, n3Entries);

    OpWrapAirSurface -> SetMaterialPropertiesTable(WrapAirProperty);
}




  // ################ LEAP ODM98 adapter (Teflon-like) ################
  // Air-Teflon Surface
    OpAirTeflonSurface = new G4OpticalSurface("Air-TeflonSurface");  
    OpAirTeflonSurface -> SetType(dielectric_metal); 
    OpAirTeflonSurface -> SetFinish(polished);
    OpAirTeflonSurface -> SetModel(glisur);

    // ODM98: 94% reflectivity guessed from datasheet for 0.75mm thickness
    G4double reflTeflon[n3Entries] = {0.94, 0.94, 0.94};   //not absorbed fraction = reflected (dielectric_metal)

    G4MaterialPropertiesTable *AirTeflonProperty = new G4MaterialPropertiesTable();
    AirTeflonProperty -> AddProperty("REFLECTIVITY", Photon3Energy, reflTeflon, n3Entries);

    OpAirTeflonSurface -> SetMaterialPropertiesTable(AirTeflonProperty);
 
  //-------- Teflon-Air Surface
    OpTeflonAirSurface = new G4OpticalSurface("Teflon-AirSurface");  
    OpTeflonAirSurface -> SetType(dielectric_metal);
    OpTeflonAirSurface -> SetFinish(polished);
    OpTeflonAirSurface -> SetModel(glisur);

    G4double reflAirTeflon[n3Entries] = {0.0, 0.0, 0.0};

    G4MaterialPropertiesTable *TeflonAirProperty = new G4MaterialPropertiesTable();
    TeflonAirProperty -> AddProperty("REFLECTIVITY", Photon3Energy, reflAirTeflon, n3Entries);
    OpTeflonAirSurface -> SetMaterialPropertiesTable(TeflonAirProperty);




    //defining the sigma_alpha for the optical surface between the scintillator and the optical coupling pad
  OpPadSciSurface = new G4OpticalSurface("Pad-SciSurface");
  OpPadSciSurface -> SetType(dielectric_dielectric);
  OpPadSciSurface -> SetModel(unified);
  OpPadSciSurface -> SetFinish(ground);
  OpPadSciSurface -> SetSigmaAlpha(sigma_alpha);
  
  /*G4MaterialPropertiesTable* PadSciSurfProp = new G4MaterialPropertiesTable();
 
  G4double PadSciRindex[n3Entries] = {1.58, 1.58, 1.58};
  G4double PadSciSpecularlobe[n3Entries] = {1.0, 1.0, 1.0};	//respect to normal of microfacet
  G4double PadSciSpecularspike[n3Entries] = {0.0, 0.0, 0.0};	//respect to average normal
  G4double PadSciBackscatter[n3Entries] = {0.0, 0.0, 0.0};

  PadSciSurfProp -> AddProperty("RINDEX",Photon3Energy,PadSciRindex,n3Entries);
  PadSciSurfProp -> AddProperty("SPECULARLOBECONSTANT",Photon3Energy,PadSciSpecularlobe,n3Entries);
  PadSciSurfProp -> AddProperty("SPECULARSPIKECONSTANT",Photon3Energy,PadSciSpecularspike,n3Entries);
  PadSciSurfProp -> AddProperty("BACKSCATTERCONSTANT",Photon3Energy,PadSciBackscatter,n3Entries);
  
  OpPadSciSurface -> SetMaterialPropertiesTable(PadSciSurfProp);*/

  OpSciPadSurface = new G4OpticalSurface("Sci-PadSurface");
  OpSciPadSurface -> SetType(dielectric_dielectric);
  OpSciPadSurface -> SetModel(unified);
  OpSciPadSurface -> SetFinish(ground);
  OpSciPadSurface -> SetSigmaAlpha(sigma_alpha);
  
  /*G4MaterialPropertiesTable* SciPadSurfProp = new G4MaterialPropertiesTable();
 
  G4double SciPadRindex[n3Entries] = {1.0, 1.0, 1.0};//{1.5324, 1.5324, 1.5324};//
  G4double SciPadSpecularlobe[n3Entries] = {1.0, 1.0, 1.0};	//respect to normal of microfacet
  G4double SciPadSpecularspike[n3Entries] = {0.0, 0.0, 0.0};	//respect to average normal
  G4double SciPadBackscatter[n3Entries] = {0.0, 0.0, 0.0};

  SciPadSurfProp -> AddProperty("RINDEX",Photon3Energy,SciPadRindex,n3Entries);
  SciPadSurfProp -> AddProperty("SPECULARLOBECONSTANT",Photon3Energy,SciPadSpecularlobe,n3Entries);
  SciPadSurfProp -> AddProperty("SPECULARSPIKECONSTANT",Photon3Energy,SciPadSpecularspike,n3Entries);
  SciPadSurfProp -> AddProperty("BACKSCATTERCONSTANT",Photon3Energy,SciPadBackscatter,n3Entries);
   
  OpSciPadSurface -> SetMaterialPropertiesTable(SciPadSurfProp);*/

//--------  Scintillator -> Air Surface: Unified model -- Option 2
  OpSciAirSurface = new G4OpticalSurface("Sci-AirSurface");
  
  /*
  G4double sigma_alpha = 0.0036*deg;  	//Value of "good" side of scintillator
//   G4double sigma_alpha = 0.0275*deg;  	//Value of "bad" side of scintillator
//   G4double sigma_alpha = 0.5275*deg;  	//Value of "bad" side of scintillator
  */
  
  OpSciAirSurface -> SetType(dielectric_dielectric);
  OpSciAirSurface -> SetModel(unified);
  OpSciAirSurface -> SetFinish(ground);
  OpSciAirSurface -> SetSigmaAlpha(sigma_alpha);

  //G4double ScintAirEnesurf[n3Entries] = {2.038*eV, 4.144*eV}; 
  G4double ScintAirSpecularlobe[n3Entries] = {1.0, 1.0, 1.0};	//respect to normal of microfacet
  G4double ScintAirSpecularspike[n3Entries] = {0.0, 0.0, 0.0};	//respect to average normal
  G4double ScintAirBackscatter[n3Entries] = {0.0, 0.0, 0.0};
  G4double ScintAirRindex[n3Entries] = {1.0, 1.0, 1.0};

  G4MaterialPropertiesTable* SciAirSurfProp = new G4MaterialPropertiesTable();

  SciAirSurfProp -> AddProperty("RINDEX",Photon3Energy,ScintAirRindex,n3Entries);
  SciAirSurfProp -> AddProperty("SPECULARLOBECONSTANT",Photon3Energy,ScintAirSpecularlobe,n3Entries);
  SciAirSurfProp -> AddProperty("SPECULARSPIKECONSTANT",Photon3Energy,ScintAirSpecularspike,n3Entries);
  SciAirSurfProp -> AddProperty("BACKSCATTERCONSTANT",Photon3Energy,ScintAirBackscatter,n3Entries);

  OpSciAirSurface -> SetMaterialPropertiesTable(SciAirSurfProp);
  
  
  
  //-------- Air -> Scintillator Surface (same as Scint -> Air, but in inverse order)
  OpAirSciSurface = new G4OpticalSurface("Air-SciSurface");

  G4double SciRefrInd2[n3Entries] = { 1.58, 1.58, 1.58};
  
  OpAirSciSurface -> SetType(dielectric_dielectric);
  OpAirSciSurface -> SetModel(unified);
  OpAirSciSurface -> SetFinish(ground);
  OpAirSciSurface -> SetSigmaAlpha(sigma_alpha);
     
  G4MaterialPropertiesTable* AirSciSurfProp = new G4MaterialPropertiesTable();

  AirSciSurfProp -> AddProperty("RINDEX",Photon3Energy,SciRefrInd2,n3Entries);  
  AirSciSurfProp -> AddProperty("SPECULARLOBECONSTANT",Photon3Energy,ScintAirSpecularlobe,n3Entries);
  AirSciSurfProp -> AddProperty("SPECULARSPIKECONSTANT",Photon3Energy,ScintAirSpecularspike,n3Entries);
  AirSciSurfProp -> AddProperty("BACKSCATTERCONSTANT",Photon3Energy,ScintAirBackscatter,n3Entries);

  OpAirSciSurface -> SetMaterialPropertiesTable(AirSciSurfProp);
  
    materials_defined_ = true;
}

void POLAR2DetectorConstruction::ConstructModule_() {
    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    // build sensitive detector: GasLog_
    if (!materials_defined_) return;
    // define ModuleBox

    /*

    // old POLAR truncated scintillators
    G4double top_taper = 5.0;     
    G4double SBzPlane[4] = {0.0-Scint_length/2.*mm, Trunc_height-Scint_length/2.*mm, (Trunc_height + 115.0-Scint_length/2.)*mm, (Trunc_height + 115.0 + 5.0-Scint_length/2.)*mm};
    G4double SBrOuter[4] = {Scint_end*mm / 2.0, Scint_thick / 2.0,Scint_thick / 2.0, top_taper / 2.0};
    G4double SBrInner[4] = {0.0, 0.0, 0.0, 0.0};
    G4Polyhedra* ScintBox = new G4Polyhedra("ScintillatorBox", 45.0*deg, 360.0*deg, 4, 4, SBzPlane, SBrInner, SBrOuter);
    */

      // scintillator bars
    G4Box* ScintBox = new G4Box("ScintillatorBox", Scint_thick/2., Scint_thick/2., Scint_length/2.);  
    ScintLog_ = new G4LogicalVolume(ScintBox, theScint, "Scint");
    G4VisAttributes* VisAtt_ScintLog = new G4VisAttributes(true, G4Color(0., 0.67, 0.87, 0.4));
    VisAtt_ScintLog->SetForceSolid(true);
    ScintLog_->SetVisAttributes(VisAtt_ScintLog);

      // scintillator block to test the scintillation efficiency
    G4Box* ScintBlock = new G4Box("ScintillatorBlock", 500/2., 500/2., 500/2.);  
    ScintBlockLog_ = new G4LogicalVolume(ScintBlock, theScint, "Scint");
    ScintBlockLog_->SetVisAttributes(VisAtt_ScintLog);


      // LEAP Polarimeter Detector Element
    G4Box* LEAPScintBox = new G4Box("LEAPScintillatorBox", LEAPScint_thick/2., LEAPScint_thick/2., LEAPScint_length/2.);  
    if (LEAPScint_conical) {
      G4Cons* LEAPScintCone = new G4Cons("LEAPScintillatorCone", 0, 0, 0, 2*(LEAPScint_length+2*sqrt(pow(2*LEAP_PMTEntranceRadius,2)-pow(LEAP_PMTEntranceRadius,2)))/sqrt(3)/2., (LEAPScint_length+2*sqrt(pow(2*LEAP_PMTEntranceRadius,2)-pow(LEAP_PMTEntranceRadius,2)))/2., 0, 2*pi);  
      G4IntersectionSolid* LEAPScintillator = new G4IntersectionSolid("LEAPScintillator", LEAPScintBox, LEAPScintCone);
      LEAPScintLog_ = new G4LogicalVolume(LEAPScintillator, theScint, "LEAPScint");
    }
    else {
      LEAPScintLog_ = new G4LogicalVolume(LEAPScintBox, theScint, "LEAPScint");
    }
    LEAPScintLog_->SetVisAttributes(VisAtt_ScintLog);

    // test cylindrical scintillator of the size of the PMT window for LEAP
    G4Tubs* LEAPCylindricalScintBox = new G4Tubs("LEAPScintillatorBox", 0, LEAP_PMTEntranceRadius*mm, LEAPScint_length/2.*mm, 0, 2*pi);  
    LEAPCylindricalScintLog_ = new G4LogicalVolume(LEAPCylindricalScintBox, theScint, "LEAPScint");
    LEAPCylindricalScintLog_->SetVisAttributes(VisAtt_ScintLog);


    // test for Federico's neutrino detector

    G4Box* NeutrinoScintBar = new G4Box("NeutrinoScintBar", 3/2., 3/2., 120/2.);  
    NeutrinoScintLog_ = new G4LogicalVolume(NeutrinoScintBar, theScint, "NeutrinoScint");
    NeutrinoScintLog_->SetVisAttributes(VisAtt_ScintLog);


      // Grid long parts
    G4Box* s_grid_long_thin = new G4Box("grid", 0.2/2.*mm, 50.7/2.*mm, Grid_height/2.);
    l_grid_long_thin = new G4LogicalVolume(s_grid_long_thin, Accura25, "grid",0,0,0); 
    G4Box* s_grid_long_thick = new G4Box("grid", 0.5/2.*mm, 50.7/2.*mm, Grid_height/2.);
    l_grid_long_thick = new G4LogicalVolume(s_grid_long_thick, Accura25, "grid",0,0,0); 
      // Grid short parts
    G4Box* s_grid_short_thin = new G4Box("grid", 6.0/2.*mm, 0.2/2.*mm, Grid_height/2.);
    l_grid_short_thin = new G4LogicalVolume(s_grid_short_thin, Accura25, "grid",0,0,0);    
    G4Box* s_grid_short_thick = new G4Box("grid", 6.0/2.*mm, 0.5/2.*mm, Grid_height/2.);
    l_grid_short_thick = new G4LogicalVolume(s_grid_short_thick, Accura25, "grid",0,0,0); 
      // Middle cross on the top side
    G4Box* s_cross_long = new G4Box("cross", 0.5/2.*mm, 45/2.*mm, Grid_height/2.);
    l_cross_long = new G4LogicalVolume(s_cross_long, Accura25, "cross",0,0,0);
    G4Box* s_cross_short = new G4Box("cross", (45-0.5)/2/2.*mm, 0.5/2.*mm, Grid_height/2.);
    l_cross_short = new G4LogicalVolume(s_cross_short, Accura25, "cross",0,0,0);

    G4Box* s_grid_plate = new G4Box("grid", 50.3/2.*mm, 50.3/2.*mm, 0.2/2.);  
    grid_plate_test = new G4LogicalVolume(s_grid_plate, Accura25, "grid",0,0,0);


    G4VisAttributes* VisAtt_GridLog = new G4VisAttributes(true, G4Color(0.05, 0.05, 0.05, 0.5));
    VisAtt_GridLog->SetForceSolid(true);
    l_grid_short_thin->SetVisAttributes(VisAtt_GridLog);
    l_grid_short_thick->SetVisAttributes(VisAtt_GridLog);
    l_grid_long_thin->SetVisAttributes(VisAtt_GridLog);
    l_grid_long_thick->SetVisAttributes(VisAtt_GridLog);
    l_cross_long->SetVisAttributes(VisAtt_GridLog);  
    l_cross_short->SetVisAttributes(VisAtt_GridLog);  

    grid_plate_test->SetVisAttributes(VisAtt_GridLog);  

    G4double wrapping_overlap = 0;
    if (vikuiti_is_metal) {wrapping_overlap = 0.0001;}
      // Wrapping
    G4Box* vikuit_x = new G4Box("vikuit", Vikuiti_thickness+0.05+Scint_thick/2.*mm, wrapping_overlap+Vikuiti_thickness/2.*mm, (Scint_length-Grid_height)/2.);
    vikuit_x_log = new G4LogicalVolume(vikuit_x, VikuitiESR_, "vikuit",0,0,0);
    G4Box* vikuit_y = new G4Box("vikuit", wrapping_overlap+Vikuiti_thickness/2.*mm, 0.05+Scint_thick/2.*mm, (Scint_length-Grid_height)/2.);
    vikuit_y_log = new G4LogicalVolume(vikuit_y, VikuitiESR_, "vikuit",0,0,0);
    G4Box* vikuit_z = new G4Box("vikuit", 50.3/2.*mm, 50.3/2.*mm, Vikuiti_thickness/2.);  
    vikuit_z_log = new G4LogicalVolume(vikuit_z, VikuitiESR_, "vikuit",0,0,0);

    G4Box* PerfRefl_x = new G4Box("PerfRefl", Vikuiti_thickness+0.05+Scint_thick/2.*mm, 0.01+Vikuiti_thickness/2.*mm, Scint_length/2.);    // -Grid_height on z if grid
    PerfRefl_x_log = new G4LogicalVolume(PerfRefl_x, VikuitiESR_, "PerfRefl",0,0,0);
    G4Box* PerfRefl_y = new G4Box("PerfRefl", 0.01+Vikuiti_thickness/2.*mm, 0.05+Scint_thick/2.*mm, Scint_length/2.);  // -Grid_height on z if grid
    PerfRefl_y_log = new G4LogicalVolume(PerfRefl_y, VikuitiESR_, "PerfRefl",0,0,0);
    G4Box* PerfRefl_z = new G4Box("PerfRefl", 50.3/2.*mm, 50.3/2.*mm, Vikuiti_thickness/2.);  
    PerfRefl_z_log = new G4LogicalVolume(PerfRefl_z, VikuitiESR_, "PerfRefl",0,0,0);

    // LEAP Vikuiti
    G4double LEAPvikuiti_halflength;
    if (LEAPScint_conical){
      LEAPvikuiti_halflength = 0.05+LEAPScint_length/2.-(LEAPTefloncone_thick-LEAPTeflonsquare_thick)/2.;
    } else{
      LEAPvikuiti_halflength = 0.05+LEAPScint_length/2.;
    }
    G4Box* LEAPvikuit_x = new G4Box("vikuit", Vikuiti_thickness+0.05+LEAPScint_thick/2.*mm, wrapping_overlap+Vikuiti_thickness/2.*mm, LEAPvikuiti_halflength);
    LEAPvikuit_x_log = new G4LogicalVolume(LEAPvikuit_x, VikuitiESR_, "vikuit",0,0,0);
    G4Box* LEAPvikuit_y = new G4Box("vikuit", wrapping_overlap+Vikuiti_thickness/2.*mm, 0.05+LEAPScint_thick/2.*mm, LEAPvikuiti_halflength);
    LEAPvikuit_y_log = new G4LogicalVolume(LEAPvikuit_y, VikuitiESR_, "vikuit",0,0,0);
    G4Box* LEAPvikuit_z = new G4Box("vikuit", Vikuiti_thickness+0.05+LEAPScint_thick/2.*mm, Vikuiti_thickness+0.05+LEAPScint_thick/2.*mm, Vikuiti_thickness/2.);  
    LEAPvikuit_z_log = new G4LogicalVolume(LEAPvikuit_z, VikuitiESR_, "vikuit",0,0,0);

    G4Tubs* LEAPCylindricalvikuit_xy = new G4Tubs("vikuit", LEAP_PMTEntranceRadius+0.05*mm, wrapping_overlap+Vikuiti_thickness+LEAP_PMTEntranceRadius+0.05*mm, (0.05+LEAPScint_length)/2.*mm, 0, 2*pi);  
    LEAPCylindricalvikuit_xy_log = new G4LogicalVolume(LEAPCylindricalvikuit_xy, VikuitiESR_, "vikuit",0,0,0);
    G4Tubs* LEAPCylindricalvikuit_z = new G4Tubs("vikuit", 0, Vikuiti_thickness+0.05+LEAP_PMTEntranceRadius*mm, Vikuiti_thickness/2.*mm, 0, 2*pi);  
    LEAPCylindricalvikuit_z_log = new G4LogicalVolume(LEAPCylindricalvikuit_z, VikuitiESR_, "vikuit",0,0,0);

    G4VisAttributes* VisAtt_Vikuiti = new G4VisAttributes(true, G4Color(0.4, 0.4, 0.5, 0.5));
    VisAtt_Vikuiti->SetForceSolid(true);
    vikuit_x_log->SetVisAttributes(VisAtt_Vikuiti);  
    vikuit_y_log->SetVisAttributes(VisAtt_Vikuiti);  
    vikuit_z_log->SetVisAttributes(VisAtt_Vikuiti);  

    PerfRefl_x_log->SetVisAttributes(VisAtt_Vikuiti);  
    PerfRefl_y_log->SetVisAttributes(VisAtt_Vikuiti);  
    PerfRefl_z_log->SetVisAttributes(VisAtt_Vikuiti);  

    LEAPvikuit_x_log->SetVisAttributes(VisAtt_Vikuiti);  
    LEAPvikuit_y_log->SetVisAttributes(VisAtt_Vikuiti);  
    LEAPvikuit_z_log->SetVisAttributes(VisAtt_Vikuiti);  
    LEAPCylindricalvikuit_xy_log->SetVisAttributes(VisAtt_Vikuiti);  
    LEAPCylindricalvikuit_z_log->SetVisAttributes(VisAtt_Vikuiti);  


/*    // old verison with claryl + option for coating the plastic grid
    G4Box* claryl_x;
    G4Box* claryl_y;
    G4Box* claryl_z;
    if (use_coating){
      claryl_x = new G4Box("claryl", alutick/2., 5.94/2.*mm, Scint_length/2.-NoWrap_height*mm);  
      claryl_y = new G4Box("claryl", 5.96/2.*mm, alutick/2., Scint_length/2.-NoWrap_height*mm);  
      claryl_z = new G4Box("claryl", 18.8/2.*mm,18.8/2.*mm,alutick/2.);
    } else{
      claryl_x = new G4Box("claryl", alutick/2., 5.94/2.*mm, Scint_length/2.*mm);  
      claryl_y = new G4Box("claryl", 5.96/2.*mm, alutick/2., Scint_length/2.*mm);  
      claryl_z = new G4Box("claryl", 18.8/2.*mm,18.8/2.*mm,alutick/2.);
    }  

    claryl_x_log = new G4LogicalVolume(claryl_x, claryl, "claryl",0,0,0);
    claryl_y_log = new G4LogicalVolume(claryl_y, claryl, "claryl",0,0,0);
    claryl_z_log = new G4LogicalVolume(claryl_z, claryl, "claryl",0,0,0);

    G4Box* air_x;
    G4Box* air_inner_y;
    G4Box* air_outer_y;
    if (use_coating){
      air_x = new G4Box("air", Air_thick_wrap/2., 5.94/2.*mm, Scint_length/2.-NoWrap_height*mm);  
      air_inner_y = new G4Box("air", 5.94/2.*mm, Air_thick_wrap/2., Scint_length/2.-NoWrap_height*mm);  
      air_outer_y = new G4Box("air", 5.96/2.*mm, Air_thick_wrap/2., Scint_length/2.-NoWrap_height*mm);  
    } else{
      air_x = new G4Box("air", Air_thick_wrap/2., 5.94/2.*mm, Scint_length/2.*mm);  
      air_inner_y = new G4Box("air", 5.94/2.*mm, Air_thick_wrap/2., Scint_length/2.*mm);  
      air_outer_y = new G4Box("air", 5.96/2.*mm, Air_thick_wrap/2., Scint_length/2.*mm);  
    }



    G4VisAttributes* VisAtt_ClarylLog = new G4VisAttributes(true, G4Color(0.3, 0.3, 0.3, 0.9));
    VisAtt_ClarylLog->SetForceSolid(true);
    claryl_x_log->SetVisAttributes(VisAtt_ClarylLog);
    claryl_y_log->SetVisAttributes(VisAtt_ClarylLog);
    claryl_z_log->SetVisAttributes(VisAtt_ClarylLog);

    // Scintillator end coating for cross talk through grid
    if (use_coating){
      G4Box* reflcoating_x = new G4Box("reflcoating", coat_thick/2.*mm, 5.9/2.*mm, Coating_height/2.*mm);  
      G4Box* reflcoating_y = new G4Box("reflcoating", 5.9/2.+coat_thick*mm, coat_thick/2.*mm, Coating_height/2.*mm);  

      reflcoating_x_log = new G4LogicalVolume(reflcoating_x, reflcoating, "reflcoating",0,0,0);
      reflcoating_y_log = new G4LogicalVolume(reflcoating_y, reflcoating, "reflcoating",0,0,0);

      G4VisAttributes* VisAtt_reflcoatingLog = new G4VisAttributes(true, G4Color(0.35, 0.35, 0.35, 0.9));
      VisAtt_reflcoatingLog->SetForceSolid(true);
      reflcoating_x_log->SetVisAttributes(VisAtt_reflcoatingLog);
      reflcoating_y_log->SetVisAttributes(VisAtt_reflcoatingLog);
    }
    */


    G4Box* air_inner_x = new G4Box("air", 6.0/2.*mm, 0.05/2.*mm, Scint_length/2.*mm);
    air_inner_x_log = new G4LogicalVolume(air_inner_x, Environment_, "air",0,0,0);
    G4Box* air_inner_y = new G4Box("air", 0.05/2.*mm, 6.0/2.*mm, Scint_length/2.*mm);
    air_inner_y_log = new G4LogicalVolume(air_inner_y, Environment_, "air",0,0,0);
    G4Box* air_outer_x = new G4Box("air", 6.13/2.*mm, 0.02/2.*mm, Scint_length/2.*mm);
    air_outer_x_log = new G4LogicalVolume(air_outer_x, Environment_, "air",0,0,0);
    G4Box* air_outer_y = new G4Box("air", 0.02/2.*mm, 6.0/2.*mm, Scint_length/2.*mm);
    air_outer_y_log = new G4LogicalVolume(air_outer_y, Environment_, "air",0,0,0);

    G4Box* air_z = new G4Box("air", 50.3/2.*mm, 50.3/2.*mm, 0.05/2.*mm);
    air_z_log = new G4LogicalVolume(air_z, Environment_, "air",0,0,0);

    // LEAP air layers

    G4double LEAPair_halflength;
    if (LEAPScint_conical){
      LEAPair_halflength = (0.05+LEAPScint_length-LEAPTefloncone_thick+LEAPTeflonsquare_thick)/2.*mm;
    } else{
      LEAPair_halflength = 0.05+LEAPScint_length/2.*mm;
    }

    G4Box* LEAPair_x = new G4Box("air", (LEAPScint_thick+0.1)/2.*mm, 0.05/2.*mm, LEAPair_halflength);
    LEAPair_x_log = new G4LogicalVolume(LEAPair_x, Environment_, "air",0,0,0);
    G4Box* LEAPair_y = new G4Box("air", 0.05/2.*mm, (LEAPScint_thick+0.1)/2.*mm, LEAPair_halflength);
    LEAPair_y_log = new G4LogicalVolume(LEAPair_y, Environment_, "air",0,0,0);
    G4Box* LEAPair_z = new G4Box("air", (LEAPScint_thick+0.1)/2.*mm, (LEAPScint_thick+0.1)/2.*mm, 0.05/2.*mm);
    LEAPair_z_log = new G4LogicalVolume(LEAPair_z, Environment_, "air",0,0,0);
    LEAPair_z_log->SetVisAttributes(VisAtt_Vikuiti);  

    G4Tubs* LEAPCylindricalair_xy = new G4Tubs("vikuit", LEAP_PMTEntranceRadius*mm, LEAP_PMTEntranceRadius+0.05*mm, (0.05+LEAPScint_length)/2.*mm, 0, 2*pi);  
    LEAPCylindricalair_xy_log = new G4LogicalVolume(LEAPCylindricalair_xy, Environment_, "air",0,0,0);
    G4Tubs* LEAPCylindricalair_z = new G4Tubs("vikuit", 0, 0.05+LEAP_PMTEntranceRadius*mm, 0.05/2.*mm, 0, 2*pi);  
    LEAPCylindricalair_z_log = new G4LogicalVolume(LEAPCylindricalair_z, Environment_, "air",0,0,0);


    G4Tubs* TeflonOptPadapperture;
    if (LEAPScint_conical){
      TeflonOptPadapperture = new G4Tubs("TeflonOptPadapperture", 0, LEAP_PMTEntranceRadius*mm, LEAPTefloncone_thick/2.*mm, 0, 2*pi);
      //G4Cons* LEAPairlayer_teflon = new G4Cons("air", LEAP_PMTEntranceRadius, LEAP_PMTEntranceRadius+0.05/sqrt(3)*mm, LEAP_PMTEntranceRadius+(LEAPTefloncone_thick-LEAPTeflonsquare_thick-0.05)/sqrt(3)*mm, LEAP_PMTEntranceRadius+(LEAPTefloncone_thick-LEAPTeflonsquare_thick)/sqrt(3)*mm, (LEAPTefloncone_thick-LEAPTeflonsquare_thick-0.05)/2.*mm, 0, 2*pi);  
      G4Box* LEAPairlayer_teflon_box = new G4Box("air", LEAPScint_thick/2., LEAPScint_thick/2., sqrt(3)*(LEAPScint_thick/sqrt(2)-4)/2.*mm);  
      G4Cons* LEAPairlayer_teflon_cone = new G4Cons("air", LEAP_PMTEntranceRadius, LEAP_PMTEntranceRadius+0.05/sqrt(3)*mm, 17/sqrt(2)*mm, 17/sqrt(2)+0.05/sqrt(3)*mm, sqrt(3)*(LEAPScint_thick/sqrt(2)-4)/2.*mm, 0, 2*pi);  
      G4IntersectionSolid* LEAPairlayer_teflon = new G4IntersectionSolid("airlayer", LEAPairlayer_teflon_cone, LEAPairlayer_teflon_box);
      LEAPair_teflonconical_log = new G4LogicalVolume(LEAPairlayer_teflon, Environment_, "air",0,0,0); 

      G4Box* box_topappertureteflon = new G4Box("air", (LEAPScint_thick+0.1)/2.*mm, (LEAPScint_thick+0.1)/2.*mm, 0.05/2.*mm);
      G4Tubs* cone_topappertureteflon = new G4Tubs("TeflonConeTopapperture", 0, LEAP_PMTEntranceRadius+(LEAPTefloncone_thick-LEAPTeflonsquare_thick)/sqrt(3)*mm, 0.05/2.*mm, 0, 2*pi);
      G4SubtractionSolid* LEAPair_teflonsquare = new G4SubtractionSolid("airlayer", box_topappertureteflon, cone_topappertureteflon);
      LEAPair_teflonsquare_log = new G4LogicalVolume(LEAPair_teflonsquare, Environment_, "air",0,0,0);
    } else{
      TeflonOptPadapperture = new G4Tubs("TeflonOptPadapperture", 0, LEAP_PMTEntranceRadius*mm, LEAPTeflonsquare_thick/2.*mm, 0, 2*pi);
      G4Box* LEAPairlayer_teflon = new G4Box("air", LEAPScint_thick/2., LEAPScint_thick/2., 0.05/2.*mm);
      G4SubtractionSolid* LEAPair_teflonsquare = new G4SubtractionSolid("airlayer", LEAPairlayer_teflon, TeflonOptPadapperture);
      LEAPair_teflonsquare_log = new G4LogicalVolume(LEAPair_teflonsquare, Environment_, "air",0,0,0);
    }


    /*
      // to check the air layer placements
    G4VisAttributes* VisAtt_air = new G4VisAttributes(true, G4Color(0.9, 0.8, 0.5, 0.1));
    VisAtt_air->SetForceSolid(true);
    air_inner_x_log->SetVisAttributes(VisAtt_air);  
    air_inner_y_log->SetVisAttributes(VisAtt_air);
    air_outer_x_log->SetVisAttributes(VisAtt_air);  
    air_outer_y_log->SetVisAttributes(VisAtt_air);  
    air_z_log->SetVisAttributes(VisAtt_air);  
    */

    
      // Debugging purposes
    G4Box* air_test = new G4Box("air", 50.3/2.*mm,50.3/2.*mm,2/2.*mm);
    G4Box* airboxtest = new G4Box("air", 50/2.*mm, 50/2.*mm, 50/2.*mm);
    air_test_log = new G4LogicalVolume(air_test, Environment_, "air",0,0,0);
    airbbox = new G4LogicalVolume(airboxtest, Environment_, "air",0,0,0);

    G4VisAttributes* VisAtt_airtest = new G4VisAttributes(true, G4Color(0.9, 0.8, 0.5, 0.1));
    VisAtt_airtest->SetForceSolid(true);
    air_test_log->SetVisAttributes(VisAtt_airtest);  
    airbbox->SetVisAttributes(VisAtt_airtest);  
    


    // Optical pad (optional)
    G4Box* Optpad_Box = new G4Box("OpticalPad", 50.1/2.*mm, 50.1/2.*mm, Pad_height/2.*mm);
    opt_pad_log = new G4LogicalVolume(Optpad_Box, dc93_500, "OpticalPad",0,0,0);    
    //opt_pad_log = new G4LogicalVolume(Optpad_Box, Environment_, "OpticalPad",0,0,0);  //air instead of pad  
    G4VisAttributes* VisAtt_opt_pad = new G4VisAttributes(true, G4Color(0.996, 0.827, 0.027, 0.6));
    VisAtt_opt_pad->SetForceSolid(true);
    opt_pad_log->SetVisAttributes(VisAtt_opt_pad);   

    // LEAP optical pad - thickness of 0.8mm, 1mm with 20% compression
    G4Tubs* LEAPOptpad_Box = new G4Tubs("LEAPOpticalPad", 0, LEAP_PMTEntranceRadius*mm, LEAP_Optpad_thick/2.*mm, 0, 2*pi);  
    LEAPopt_pad_log = new G4LogicalVolume(LEAPOptpad_Box, dc93_500, "OpticalPad",0,0,0);    
    LEAPopt_pad_log->SetVisAttributes(VisAtt_opt_pad);   


    // SiPM coating
    G4Box* EpoxyResin_Box = new G4Box("EpoxyResin", 25/2.*mm, 25/2.*mm, EpoxyResin_height/2.);
    EpoxyResin_log = new G4LogicalVolume(EpoxyResin_Box, EpoxyResin, "EpoxyResin",0,0,0);    
    G4VisAttributes* VisAtt_EpoxyResin = new G4VisAttributes(true, G4Color(0.56, 0.59, 0.47, 0.8));
    VisAtt_EpoxyResin->SetForceSolid(true);
    EpoxyResin_log->SetVisAttributes(VisAtt_EpoxyResin);    

    // SiPM
    G4Box* SiPMBox = new G4Box("SiPMBox", 6.0/2.*mm, 6.0/2.*mm, SiPMsilicon_height/2.*mm);
    SiPMLog_ = new G4LogicalVolume(SiPMBox, FR4_, "SiPM",0,0,0);
    G4VisAttributes* VisAtt_SiPMLog = new G4VisAttributes(true, G4Color(0.56, 0.59, 0.47, 1.0));
    VisAtt_SiPMLog->SetForceSolid(true);
    SiPMLog_->SetVisAttributes(VisAtt_SiPMLog);
   
    // LEAP PMT
    G4Tubs* LEAPPMTBox = new G4Tubs("LEAPPMT", 0, LEAP_PMTRadius*mm, LEAP_PMT_thick/2.*mm, 0, 2*pi);  
    LEAPPMTLog_ = new G4LogicalVolume(LEAPPMTBox, FR4_, "LEAPPMT",0,0,0);
    LEAPPMTLog_->SetVisAttributes(VisAtt_SiPMLog);

    // LEAP Teflon adapter

    G4Box* TeflonBox;
    if (LEAPScint_conical){
      TeflonBox = new G4Box("TeflonBox", (0.3+LEAPScint_thick)/2.*mm, (0.3+LEAPScint_thick)/2.*mm, LEAPTefloncone_thick/2.*mm); //300um bigger to enclose the scint with the vikuiti
    }else{
      TeflonBox = new G4Box("TeflonBox", (0.3+LEAPScint_thick)/2.*mm, (0.3+LEAPScint_thick)/2.*mm, LEAPTeflonsquare_thick/2.*mm);  // vikuiti goes around the teflon piece, to be modified
    }

    // hole for opt pad
    G4SubtractionSolid* Teflon_subOptPad = new G4SubtractionSolid("TeflonAdapter", TeflonBox, TeflonOptPadapperture);

    // hole for PMT
    G4RotationMatrix rotm = G4RotationMatrix();
    G4Transform3D trPMTapperture;
    if (LEAPScint_conical){
      trPMTapperture = G4Transform3D(rotm, G4ThreeVector(0.,0.,-(LEAPTefloncone_thick-(LEAPTeflonsquare_thick-LEAP_Optpad_apperture_thick))/2.));
    } else{
      trPMTapperture = G4Transform3D(rotm, G4ThreeVector(0.,0.,-LEAP_Optpad_apperture_thick/2.));
    }
    G4Tubs* TeflonPMTapperture = new G4Tubs("TeflonPMTapperture", 0, LEAP_PMTRadius*mm, (LEAPTeflonsquare_thick-LEAP_Optpad_apperture_thick)/2.*mm, 0, 2*pi);  
    G4SubtractionSolid* TeflonAdapter_square = new G4SubtractionSolid("TeflonAdapter", Teflon_subOptPad, TeflonPMTapperture, trPMTapperture);

    // hole for cone if conical is true -> here
    G4Transform3D trConeapperture = G4Transform3D(rotm, G4ThreeVector(0.,0.,LEAPTeflonsquare_thick/2.));
    G4Cons* TeflonConeapperture = new G4Cons("TeflonConeapperture", 0, LEAP_PMTEntranceRadius*mm, 0, LEAP_PMTEntranceRadius+(LEAPTefloncone_thick-LEAPTeflonsquare_thick)/sqrt(3)*mm, (LEAPTefloncone_thick-LEAPTeflonsquare_thick)/2.*mm, 0, 2*pi);  
    G4SubtractionSolid* TeflonAdapter = new G4SubtractionSolid("TeflonAdapter", TeflonAdapter_square, TeflonConeapperture, trConeapperture);


    //G4SubtractionSolid* TeflonAdapter = new G4SubtractionSolid("TeflonAdapter", TeflonBox, TeflonApperture);
    if (LEAPScint_conical){
      TeflonLog_ = new G4LogicalVolume(TeflonAdapter, VikuitiESR_, "TeflonAdapter",0,0,0);  // ODM98 material is actually not Teflon, density of 1.5g/cm3, we put same material as Vikuiti here
    } else{
      TeflonLog_ = new G4LogicalVolume(TeflonAdapter_square, VikuitiESR_, "TeflonAdapter",0,0,0);
    }

    G4VisAttributes* VisAtt_TeflonLog = new G4VisAttributes(true, G4Color(0.7, 0.7, 0.7, 0.4));
    VisAtt_TeflonLog->SetForceSolid(true);
    TeflonLog_->SetVisAttributes(VisAtt_TeflonLog);


    /*G4Tubs* TeflonPMTapperture = new G4Tubs("TeflonPMTapperture", 0, LEAP_PMTRadius*mm, (LEAPTeflonsquare_thick-LEAP_Optpad_apperture_thick)/2.*mm, 0, 2*pi);  
    G4Tubs* TeflonOptPadapperture = new G4Tubs("TeflonOptPadapperture", 0, LEAP_PMTEntranceRadius*mm, LEAP_Optpad_apperture_thick*mm, 0, 2*pi);

    G4RotationMatrix rotm = G4RotationMatrix();
    G4Transform3D trPMTapperture = G4Transform3D(rotm, G4ThreeVector(0.,0.,-LEAP_Optpad_apperture_thick/2.));
    G4Transform3D trOptPadapperture = G4Transform3D(rotm, G4ThreeVector(0.,0.,(LEAPTeflonsquare_thick-LEAP_Optpad_apperture_thick)/2.));

    G4MultiUnion* TeflonApperture = new G4MultiUnion("TeflonApperture");
    TeflonApperture->AddNode(*TeflonOptPadapperture, trOptPadapperture); // add shapes to the structure
    TeflonApperture->AddNode(*TeflonPMTapperture, trPMTapperture);
    //if (LEAPScint_conical){
      // define here the conical hole in teflon and add it to multiunion
      //G4Cons* TeflonConicalapperture = ;
      //TeflonApperture->AddNode(*TeflonConicalapperture,);
    //}
    TeflonApperture->Voxelize(); // close the structure

    G4SubtractionSolid* TeflonAdapter = new G4SubtractionSolid("TeflonAdapter", TeflonBox, TeflonApperture);
    TeflonLog_ = new G4LogicalVolume(TeflonAdapter, VikuitiESR_, "TeflonAdapter",0,0,0);  // teflon material + properties need to be defined

    G4VisAttributes* VisAtt_TeflonLog = new G4VisAttributes(true, G4Color(0.7, 0.7, 0.7, 1.0));
    VisAtt_TeflonLog->SetForceSolid(true);
    TeflonLog_->SetVisAttributes(VisAtt_TeflonLog);*/


   
}

void POLAR2DetectorConstruction::ConstructDetector_() {
    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    // define DetectorBox
    
    G4double det_hx = 500.0*mm;
    G4double det_hy = 500.0*mm;
    G4double det_hz = 500.0*mm;
    G4Box* DetectorBox = new G4Box("DetectorBox", det_hx, det_hy, det_hz);
    
    DetectorLog_ = new G4LogicalVolume(DetectorBox, Environment_, "DetectorLog");
    G4VisAttributes* VisAtt_DetectorLog = new G4VisAttributes(true, G4Colour(1.0, 0.0, 0.0, 1.0));
    VisAtt_DetectorLog->SetForceWireframe(true);
    DetectorLog_->SetVisAttributes(VisAtt_DetectorLog);

/*    // ################# from perfect detector to real detector #################
        // Defining arrays
    G4VPhysicalVolume* scintillator_bars[nScintBars];
    G4VPhysicalVolume* sipm_channels[nScintBars];
    G4VPhysicalVolume* PerfRefl_x1[nScintBars];
    G4VPhysicalVolume* PerfRefl_x2[nScintBars];
    G4VPhysicalVolume* PerfRefl_y1[nScintBars];
    G4VPhysicalVolume* PerfRefl_y2[nScintBars];
    G4VPhysicalVolume* air_inner_x1[nScintBars];
    G4VPhysicalVolume* air_inner_x2[nScintBars];
    G4VPhysicalVolume* air_inner_y1[nScintBars];
    G4VPhysicalVolume* air_inner_y2[nScintBars];
    G4VPhysicalVolume* air_outer_x1[nScintBars];
    G4VPhysicalVolume* air_outer_x2[nScintBars];
    G4VPhysicalVolume* air_outer_y1[nScintBars];
    G4VPhysicalVolume* air_outer_y2[nScintBars];
    G4VPhysicalVolume* long_grid_thick[3];
    G4VPhysicalVolume* long_grid_thin[6];
    G4VPhysicalVolume* short_grid_thick[24];
    G4VPhysicalVolume* short_grid_thin[48];

    G4VPhysicalVolume* vikuiti_x1[nScintBars];
    G4VPhysicalVolume* vikuiti_x2[nScintBars];
    G4VPhysicalVolume* vikuiti_y1[nScintBars];
    G4VPhysicalVolume* vikuiti_y2[nScintBars];

        // Scintillator bars and wrapping
    const G4double x_start = -(BarPitch * 3 + (BarPitch + 0.3) / 2);
    const G4double y_start = -(BarPitch * 3 + (BarPitch + 0.3) / 2);
    for (int iBar = 0; iBar < nScintBars; iBar++) {
        G4double bar_x = x_start + (iBar % 8) * BarPitch;
        G4double bar_y = y_start + (iBar / 8) * BarPitch;
        if (iBar%8 > 3){bar_x += 0.3;}  // the pitch in the middle cross is 6.5mm instead of 6.2mm due to the space of 0.1mm betwenn arrays (plus we have the 200um twice)
        if (iBar/8 > 3){bar_y += 0.3;}
        scintillator_bars[iBar] = new G4PVPlacement(NULL, G4ThreeVector(bar_x, bar_y, 0), ScintLog_, "Scintillator", DetectorLog_, false, iBar);
            // Vikuiti wrapping
        PerfRefl_x1[iBar] = new G4PVPlacement(NULL, G4ThreeVector(bar_x, bar_y-(Scint_thick/2.+0.05+Vikuiti_thickness/2.), Grid_height/2.*mm), PerfRefl_x_log, "PerfRefl", DetectorLog_, false, 2*iBar);
        PerfRefl_x2[iBar] = new G4PVPlacement(NULL, G4ThreeVector(bar_x, bar_y+(Scint_thick/2.+0.05+Vikuiti_thickness/2.), Grid_height/2.*mm), PerfRefl_x_log, "PerfRefl", DetectorLog_, false, 2*iBar+1);
        PerfRefl_y1[iBar] = new G4PVPlacement(NULL, G4ThreeVector(bar_x-(Scint_thick/2.+0.05+Vikuiti_thickness/2.), bar_y, Grid_height/2.*mm), PerfRefl_y_log, "PerfRefl", DetectorLog_, false, 2*iBar);
        PerfRefl_y2[iBar] = new G4PVPlacement(NULL, G4ThreeVector(bar_x+(Scint_thick/2.+0.05+Vikuiti_thickness/2.), bar_y, Grid_height/2.*mm), PerfRefl_y_log, "PerfRefl", DetectorLog_, false, 2*iBar+1);
        //vikuiti_x1[iBar] = new G4PVPlacement(NULL, G4ThreeVector(bar_x, bar_y-(Scint_thick/2.+0.05+Vikuiti_thickness/2.), Grid_height/2.*mm), vikuit_x_log, "Vikuiti", DetectorLog_, false, 2*iBar);
        //vikuiti_x2[iBar] = new G4PVPlacement(NULL, G4ThreeVector(bar_x, bar_y+(Scint_thick/2.+0.05+Vikuiti_thickness/2.), Grid_height/2.*mm), vikuit_x_log, "Vikuiti", DetectorLog_, false, 2*iBar+1);
        //vikuiti_y1[iBar] = new G4PVPlacement(NULL, G4ThreeVector(bar_x-(Scint_thick/2.+0.05+Vikuiti_thickness/2.), bar_y, Grid_height/2.*mm), vikuit_y_log, "Vikuiti", DetectorLog_, false, 2*iBar);
        //vikuiti_y2[iBar] = new G4PVPlacement(NULL, G4ThreeVector(bar_x+(Scint_thick/2.+0.05+Vikuiti_thickness/2.), bar_y, Grid_height/2.*mm), vikuit_y_log, "Vikuiti", DetectorLog_, false, 2*iBar+1);
            // Air strips
        air_inner_x1[iBar] = new G4PVPlacement(NULL, G4ThreeVector(bar_x, bar_y-(Scint_thick/2.+0.05/2.), 0), air_inner_x_log, "Air", DetectorLog_, false, 2*iBar);
        air_inner_x2[iBar] = new G4PVPlacement(NULL, G4ThreeVector(bar_x, bar_y+(Scint_thick/2.+0.05/2.), 0), air_inner_x_log, "Air", DetectorLog_, false, 2*iBar+1);
        air_inner_y1[iBar] = new G4PVPlacement(NULL, G4ThreeVector(bar_x-(Scint_thick/2.+0.05/2.), bar_y, 0), air_inner_y_log, "Air", DetectorLog_, false, 2*iBar);
        air_inner_y2[iBar] = new G4PVPlacement(NULL, G4ThreeVector(bar_x+(Scint_thick/2.+0.05/2.), bar_y, 0), air_inner_y_log, "Air", DetectorLog_, false, 2*iBar+1);
        air_outer_x1[iBar] = new G4PVPlacement(NULL, G4ThreeVector(bar_x, bar_y-(Scint_thick/2.+0.05+Vikuiti_thickness+0.02/2.), 0), air_outer_x_log, "Air", DetectorLog_, false, 2*iBar);
        air_outer_x2[iBar] = new G4PVPlacement(NULL, G4ThreeVector(bar_x, bar_y+(Scint_thick/2.+0.05+Vikuiti_thickness+0.02/2.), 0), air_outer_x_log, "Air", DetectorLog_, false, 2*iBar+1);
        air_outer_y1[iBar] = new G4PVPlacement(NULL, G4ThreeVector(bar_x-(Scint_thick/2.+0.05+Vikuiti_thickness+0.02/2.), bar_y, 0), air_outer_y_log, "Air", DetectorLog_, false, 2*iBar);
        air_outer_y2[iBar] = new G4PVPlacement(NULL, G4ThreeVector(bar_x+(Scint_thick/2.+0.05+Vikuiti_thickness+0.02/2.), bar_y, 0), air_outer_y_log, "Air", DetectorLog_, false, 2*iBar+1);
           // SiPM channels
        sipm_channels[iBar] = new G4PVPlacement(NULL, G4ThreeVector(bar_x, bar_y, -Pad_height-EpoxyResin_height-Scint_length/2.-SiPMsilicon_height/2.*mm), SiPMLog_, "SiPM", DetectorLog_, false, iBar);
    }

        // PerfRefl top layer
    G4VPhysicalVolume* vikuiti_top = new G4PVPlacement(NULL, G4ThreeVector(0, 0, 0.05+(Scint_length+Vikuiti_thickness)/2.*mm), vikuit_z_log, "Vikuiti", DetectorLog_, false, 0);
    G4VPhysicalVolume* air_top = new G4PVPlacement(NULL, G4ThreeVector(0, 0, (Scint_length+0.05)/2.*mm), air_z_log, "Air", DetectorLog_, false, 0);
    G4VPhysicalVolume* air_top2 = new G4PVPlacement(NULL, G4ThreeVector(0, 0, (Scint_length+0.05)/2.+Vikuiti_thickness+0.05*mm), air_z_log, "Air", DetectorLog_, false, 1);

        // Alignement grid
    for (int ii = 0; ii < 3; ii++) {
        long_grid_thick[ii] = new G4PVPlacement(NULL, G4ThreeVector((ii-1)*25.1*mm, 0, -Scint_length/2.+Grid_height/2.*mm), l_grid_long_thick, "Grid", DetectorLog_, false, ii);   
        long_grid_thin[2*ii] = new G4PVPlacement(NULL, G4ThreeVector((6.35+ii*BarPitch)*mm, 0, -Scint_length/2.+Grid_height/2.*mm), l_grid_long_thin, "Grid", DetectorLog_, false, 2*ii);   
        long_grid_thin[2*ii+1] = new G4PVPlacement(NULL, G4ThreeVector(-(6.35+ii*BarPitch)*mm, 0, -Scint_length/2.+Grid_height/2.*mm), l_grid_long_thin, "Grid", DetectorLog_, false, 2*ii+1);   
        for (int jj = 0; jj < 4; jj++) {
            short_grid_thin[16*ii+4*jj] = new G4PVPlacement(NULL, G4ThreeVector(-(3.25+jj*BarPitch)*mm, -(6.35+ii*BarPitch)*mm, -Scint_length/2.+Grid_height/2.*mm), l_grid_short_thin, "Grid", DetectorLog_, false, 16*ii+4*jj);
            short_grid_thin[16*ii+4*jj+1] = new G4PVPlacement(NULL, G4ThreeVector((3.25+jj*BarPitch)*mm, -(6.35+ii*BarPitch)*mm, -Scint_length/2.+Grid_height/2.*mm), l_grid_short_thin, "Grid", DetectorLog_, false, 16*ii+4*jj+1);
            short_grid_thin[16*ii+4*jj+2] = new G4PVPlacement(NULL, G4ThreeVector(-(3.25+jj*BarPitch)*mm, (6.35+ii*BarPitch)*mm, -Scint_length/2.+Grid_height/2.*mm), l_grid_short_thin, "Grid", DetectorLog_, false, 16*ii+4*jj+2);
            short_grid_thin[16*ii+4*jj+3] = new G4PVPlacement(NULL, G4ThreeVector((3.25+jj*BarPitch)*mm, (6.35+ii*BarPitch)*mm, -Scint_length/2.+Grid_height/2.*mm), l_grid_short_thin, "Grid", DetectorLog_, false, 16*ii+4*jj+3);
            short_grid_thick[8*ii+2*jj] = new G4PVPlacement(NULL, G4ThreeVector(-(3.25+jj*BarPitch)*mm, (ii-1)*25.1*mm, -Scint_length/2.+Grid_height/2.*mm), l_grid_short_thick, "Grid", DetectorLog_, false, 8*ii+2*jj);
            short_grid_thick[8*ii+2*jj+1] = new G4PVPlacement(NULL, G4ThreeVector((3.25+jj*BarPitch)*mm, (ii-1)*25.1*mm, -Scint_length/2.+Grid_height/2.*mm), l_grid_short_thick, "Grid", DetectorLog_, false, 8*ii+2*jj+1);
        }
    }

        // Alignement cross
    new G4PVPlacement(NULL, G4ThreeVector(0, 0, Scint_length/2.-Grid_height/2.*mm), l_cross_long, "Cross", DetectorLog_, false, 0);
    new G4PVPlacement(NULL, G4ThreeVector(-45/2/2.*mm, 0, Scint_length/2.-Grid_height/2.*mm), l_cross_short, "Cross", DetectorLog_, false, 0);
    new G4PVPlacement(NULL, G4ThreeVector(45/2/2.*mm, 0, Scint_length/2.-Grid_height/2.*mm), l_cross_short, "Cross", DetectorLog_, false, 1);

        // Optical pad
    G4VPhysicalVolume* optical_pad = new G4PVPlacement(NULL, G4ThreeVector(0,0, -Scint_length/2.-Pad_height/2.*mm), opt_pad_log, "OpticalPad", DetectorLog_, false, 0);
    
        // Epoxy resin SiPM
    G4VPhysicalVolume* epoxy_resin_0 = new G4PVPlacement(NULL, G4ThreeVector(25.1/2.*mm, 25.1/2*mm, -Scint_length/2.-Pad_height-EpoxyResin_height/2.*mm), EpoxyResin_log, "EpoxyResin", DetectorLog_, false, 0);
    G4VPhysicalVolume* epoxy_resin_1 = new G4PVPlacement(NULL, G4ThreeVector(25.1/2.*mm, -25.1/2*mm, -Scint_length/2.-Pad_height-EpoxyResin_height/2.*mm), EpoxyResin_log, "EpoxyResin", DetectorLog_, false, 1);
    G4VPhysicalVolume* epoxy_resin_2 = new G4PVPlacement(NULL, G4ThreeVector(-25.1/2.*mm, 25.1/2*mm, -Scint_length/2.-Pad_height-EpoxyResin_height/2.*mm), EpoxyResin_log, "EpoxyResin", DetectorLog_, false, 2);
    G4VPhysicalVolume* epoxy_resin_3 = new G4PVPlacement(NULL, G4ThreeVector(-25.1/2.*mm, -25.1/2*mm, -Scint_length/2.-Pad_height-EpoxyResin_height/2.*mm), EpoxyResin_log, "EpoxyResin", DetectorLog_, false, 3);


        // Defining interfaces
    // top PerfRefl layer
    new G4LogicalBorderSurface("Vikuiti-AirSurface", vikuiti_top, air_top, OpWrapAirSurface);
    new G4LogicalBorderSurface("Air-VikuitiSurface", air_top, vikuiti_top, OpAirWrapSurface);
    new G4LogicalBorderSurface("Vikuiti-AirSurface", vikuiti_top, air_top2, OpWrapAirSurface);
    new G4LogicalBorderSurface("Air-VikuitiSurface", air_top2, vikuiti_top, OpAirWrapSurface);
    for (int iBar = 0; iBar < nScintBars; iBar++) {
        G4int col = iBar%8;
        G4int row = iBar/8; 
        // scintillator 4 sides
        new G4LogicalBorderSurface("Sci-AirSurface", scintillator_bars[iBar], air_inner_x1[iBar], OpSciAirSurface);
        new G4LogicalBorderSurface("Air-SciSurface", air_inner_x1[iBar], scintillator_bars[iBar], OpAirSciSurface);
        new G4LogicalBorderSurface("Sci-AirSurface", scintillator_bars[iBar], air_inner_x2[iBar], OpSciAirSurface);
        new G4LogicalBorderSurface("Air-SciSurface", air_inner_x2[iBar], scintillator_bars[iBar], OpAirSciSurface);
        new G4LogicalBorderSurface("Sci-AirSurface", scintillator_bars[iBar], air_inner_y1[iBar], OpSciAirSurface);
        new G4LogicalBorderSurface("Air-SciSurface", air_inner_y1[iBar], scintillator_bars[iBar], OpAirSciSurface);
        new G4LogicalBorderSurface("Sci-AirSurface", scintillator_bars[iBar], air_inner_y2[iBar], OpSciAirSurface);
        new G4LogicalBorderSurface("Air-SciSurface", air_inner_y2[iBar], scintillator_bars[iBar], OpAirSciSurface);
        // scintillator top
        new G4LogicalBorderSurface("Sci-AirSurface", scintillator_bars[iBar], air_top, OpSciAirSurface);
        new G4LogicalBorderSurface("Air-SciSurface", air_top, scintillator_bars[iBar], OpAirSciSurface);
        // scintillator bottom (with optical pad)
        new G4LogicalBorderSurface("Sci-PadSurface", scintillator_bars[iBar], optical_pad, OpSciPadSurface);
        new G4LogicalBorderSurface("Pad-SciSurface", optical_pad, scintillator_bars[iBar], OpPadSciSurface);
        //perfect wrapping
        //new G4LogicalBorderSurface("PerfRefl-AirSurface", PerfRefl_x1[iBar], air_inner_x1[iBar], OpPerfReflAirSurface);
        new G4LogicalBorderSurface("PerfRefl-AirSurface", PerfRefl_x1[iBar], air_outer_x1[iBar], OpPerfReflAirSurface);
        new G4LogicalBorderSurface("Air-PerfReflSurface", air_inner_x1[iBar], PerfRefl_x1[iBar], OpAirPerfReflSurface);
        new G4LogicalBorderSurface("Air-PerfReflSurface", air_outer_x1[iBar], PerfRefl_x1[iBar], OpAirPerfReflSurface);
        new G4LogicalBorderSurface("PerfRefl-AirSurface", PerfRefl_x2[iBar], air_inner_x2[iBar], OpPerfReflAirSurface);
        new G4LogicalBorderSurface("PerfRefl-AirSurface", PerfRefl_x2[iBar], air_outer_x2[iBar], OpPerfReflAirSurface);
        new G4LogicalBorderSurface("Air-PerfReflSurface", air_inner_x2[iBar], PerfRefl_x2[iBar], OpAirPerfReflSurface);
        new G4LogicalBorderSurface("Air-PerfReflSurface", air_outer_x2[iBar], PerfRefl_x2[iBar], OpAirPerfReflSurface);
        new G4LogicalBorderSurface("PerfRefl-AirSurface", PerfRefl_y1[iBar], air_inner_y1[iBar], OpPerfReflAirSurface);
        new G4LogicalBorderSurface("PerfRefl-AirSurface", PerfRefl_y1[iBar], air_outer_y1[iBar], OpPerfReflAirSurface);
        new G4LogicalBorderSurface("Air-PerfReflSurface", air_inner_y1[iBar], PerfRefl_y1[iBar], OpAirPerfReflSurface);
        new G4LogicalBorderSurface("Air-PerfReflSurface", air_outer_y1[iBar], PerfRefl_y1[iBar], OpAirPerfReflSurface);
        new G4LogicalBorderSurface("PerfRefl-AirSurface", PerfRefl_y2[iBar], air_inner_y2[iBar], OpPerfReflAirSurface);
        new G4LogicalBorderSurface("PerfRefl-AirSurface", PerfRefl_y2[iBar], air_outer_y2[iBar], OpPerfReflAirSurface);
        new G4LogicalBorderSurface("Air-PerfReflSurface", air_inner_y2[iBar], PerfRefl_y2[iBar], OpAirPerfReflSurface);
        //new G4LogicalBorderSurface("Air-PerfReflSurface", air_outer_y2[iBar], PerfRefl_y2[iBar], OpAirPerfReflSurface);
        //new G4LogicalBorderSurface("Vikuiti-AirSurface", vikuiti_x1[iBar], air_inner_x1[iBar], OpPerfReflAirSurface);
        new G4LogicalBorderSurface("Vikuiti-AirSurface", vikuiti_x1[iBar], air_outer_x1[iBar], OpPerfReflAirSurface);
        new G4LogicalBorderSurface("Air-VikuitiSurface", air_inner_x1[iBar], vikuiti_x1[iBar], OpAirPerfReflSurface);
        new G4LogicalBorderSurface("Air-VikuitiSurface", air_outer_x1[iBar], vikuiti_x1[iBar], OpAirPerfReflSurface);
        new G4LogicalBorderSurface("Vikuiti-AirSurface", vikuiti_x2[iBar], air_inner_x2[iBar], OpPerfReflAirSurface);
        new G4LogicalBorderSurface("Vikuiti-AirSurface", vikuiti_x2[iBar], air_outer_x2[iBar], OpPerfReflAirSurface);
        new G4LogicalBorderSurface("Air-VikuitiSurface", air_inner_x2[iBar], vikuiti_x2[iBar], OpAirPerfReflSurface);
        new G4LogicalBorderSurface("Air-VikuitiSurface", air_outer_x2[iBar], vikuiti_x2[iBar], OpAirPerfReflSurface);
        new G4LogicalBorderSurface("Vikuiti-AirSurface", vikuiti_y1[iBar], air_inner_y1[iBar], OpPerfReflAirSurface);
        new G4LogicalBorderSurface("Vikuiti-AirSurface", vikuiti_y1[iBar], air_outer_y1[iBar], OpPerfReflAirSurface);
        new G4LogicalBorderSurface("Air-VikuitiSurface", air_inner_y1[iBar], vikuiti_y1[iBar], OpAirPerfReflSurface);
        new G4LogicalBorderSurface("Air-VikuitiSurface", air_outer_y1[iBar], vikuiti_y1[iBar], OpAirPerfReflSurface);
        new G4LogicalBorderSurface("Vikuiti-AirSurface", vikuiti_y2[iBar], air_inner_y2[iBar], OpPerfReflAirSurface);
        new G4LogicalBorderSurface("Vikuiti-AirSurface", vikuiti_y2[iBar], air_outer_y2[iBar], OpPerfReflAirSurface);
        new G4LogicalBorderSurface("Air-VikuitiSurface", air_inner_y2[iBar], vikuiti_y2[iBar], OpAirPerfReflSurface);
        //new G4LogicalBorderSurface("Air-VikuitiSurface", air_outer_y2[iBar], vikuiti_y2[iBar], OpAirPerfReflSurface);
        //vikuiti wrapping
        new G4LogicalBorderSurface("Vikuiti-AirSurface", vikuiti_x1[iBar], air_inner_x1[iBar], OpWrapAirSurface);
        new G4LogicalBorderSurface("Vikuiti-AirSurface", vikuiti_x1[iBar], air_outer_x1[iBar], OpWrapAirSurface);
        new G4LogicalBorderSurface("Air-VikuitiSurface", air_inner_x1[iBar], vikuiti_x1[iBar], OpAirWrapSurface);
        new G4LogicalBorderSurface("Air-VikuitiSurface", air_outer_x1[iBar], vikuiti_x1[iBar], OpAirWrapSurface);
        new G4LogicalBorderSurface("Vikuiti-AirSurface", vikuiti_x2[iBar], air_inner_x2[iBar], OpWrapAirSurface);
        new G4LogicalBorderSurface("Vikuiti-AirSurface", vikuiti_x2[iBar], air_outer_x2[iBar], OpWrapAirSurface);
        new G4LogicalBorderSurface("Air-VikuitiSurface", air_inner_x2[iBar], vikuiti_x2[iBar], OpAirWrapSurface);
        new G4LogicalBorderSurface("Air-VikuitiSurface", air_outer_x2[iBar], vikuiti_x2[iBar], OpAirWrapSurface);
        new G4LogicalBorderSurface("Vikuiti-AirSurface", vikuiti_y1[iBar], air_inner_y1[iBar], OpWrapAirSurface);
        new G4LogicalBorderSurface("Vikuiti-AirSurface", vikuiti_y1[iBar], air_outer_y1[iBar], OpWrapAirSurface);
        new G4LogicalBorderSurface("Air-VikuitiSurface", air_inner_y1[iBar], vikuiti_y1[iBar], OpAirWrapSurface);
        new G4LogicalBorderSurface("Air-VikuitiSurface", air_outer_y1[iBar], vikuiti_y1[iBar], OpAirWrapSurface);
        new G4LogicalBorderSurface("Vikuiti-AirSurface", vikuiti_y2[iBar], air_inner_y2[iBar], OpWrapAirSurface);
        new G4LogicalBorderSurface("Vikuiti-AirSurface", vikuiti_y2[iBar], air_outer_y2[iBar], OpWrapAirSurface);
        new G4LogicalBorderSurface("Air-VikuitiSurface", air_inner_y2[iBar], vikuiti_y2[iBar], OpAirWrapSurface);
        new G4LogicalBorderSurface("Air-VikuitiSurface", air_outer_y2[iBar], vikuiti_y2[iBar], OpAirWrapSurface);
        // sipm PDE        
        if (row<4){
          if (col<4){
            new G4LogicalBorderSurface("ResinSiPMSurface", epoxy_resin_3, sipm_channels[iBar], OpResinSiPMSurface);
            new G4LogicalBorderSurface("ResinSiPMSurface", sipm_channels[iBar], epoxy_resin_3, OpSiPMResinSurface);
          }else{
            new G4LogicalBorderSurface("ResinSiPMSurface", epoxy_resin_1, sipm_channels[iBar], OpResinSiPMSurface);
            new G4LogicalBorderSurface("ResinSiPMSurface", sipm_channels[iBar], epoxy_resin_1, OpSiPMResinSurface);
          }
        }else{
          if (col<4){
            new G4LogicalBorderSurface("ResinSiPMSurface", epoxy_resin_2, sipm_channels[iBar], OpResinSiPMSurface);
            new G4LogicalBorderSurface("ResinSiPMSurface", sipm_channels[iBar], epoxy_resin_2, OpSiPMResinSurface);
          }else{
            new G4LogicalBorderSurface("ResinSiPMSurface", epoxy_resin_0, sipm_channels[iBar], OpResinSiPMSurface);
            new G4LogicalBorderSurface("ResinSiPMSurface", sipm_channels[iBar], epoxy_resin_0, OpSiPMResinSurface);
          }
        }
        if (row<4){
          if (col<4){
            new G4LogicalBorderSurface("ResinSiPMSurface", scintillator_bars[iBar], sipm_channels[iBar], OpResinSiPMSurface);
            new G4LogicalBorderSurface("ResinSiPMSurface", sipm_channels[iBar], scintillator_bars[iBar], OpSiPMResinSurface);
          }else{
            new G4LogicalBorderSurface("ResinSiPMSurface", scintillator_bars[iBar], sipm_channels[iBar], OpResinSiPMSurface);
            new G4LogicalBorderSurface("ResinSiPMSurface", sipm_channels[iBar], scintillator_bars[iBar], OpSiPMResinSurface);
          }
        }else{
          if (col<4){
            new G4LogicalBorderSurface("ResinSiPMSurface", scintillator_bars[iBar], sipm_channels[iBar], OpResinSiPMSurface);
            new G4LogicalBorderSurface("ResinSiPMSurface", sipm_channels[iBar], scintillator_bars[iBar], OpSiPMResinSurface);
          }else{
            new G4LogicalBorderSurface("ResinSiPMSurface", scintillator_bars[iBar], sipm_channels[iBar], OpResinSiPMSurface);
            new G4LogicalBorderSurface("ResinSiPMSurface", sipm_channels[iBar], scintillator_bars[iBar], OpSiPMResinSurface);
          }
        }
        // grid
        G4int grid_index;
        if(row==0 || row==4){ // short thick grid parts
            grid_index = (6-2*(iBar-32*(row/4)))*(1-(iBar-32*(row/4))/4)+(2*(iBar-32*(row/4)-4)+1)*((iBar-32*(row/4))/4)+8*(row/4);
            new G4LogicalBorderSurface("Grid-AirSurface", short_grid_thick[grid_index], air_inner_x1[iBar], OpGridAirSurface);
            new G4LogicalBorderSurface("Air-GridSurface", air_inner_x1[iBar], short_grid_thick[grid_index], OpAirGridSurface);
        }else{ // short thin grid parts
            grid_index = (44-16*(row-1))*(1-(row/4))+(14+16*(row-5))*(row/4)-4*(iBar-8*row)*(1-((iBar-8*row)/4))-11*((iBar-8*row)/4)+4*(iBar-8*row-4)*((iBar-8*row)/4);
            new G4LogicalBorderSurface("Grid-AirSurface", short_grid_thin[grid_index], air_inner_x1[iBar], OpGridAirSurface);
            new G4LogicalBorderSurface("Air-GridSurface", air_inner_x1[iBar], short_grid_thin[grid_index], OpAirGridSurface);
        }
        if(row==3 || row==7){ // short thick grid parts
            grid_index = (14-2*(iBar-32*(row/7)-24))*(1-((iBar-32*(row/7)-24)/4))+(2*(iBar-32*(row/7)-28)+9)*((iBar-32*(row/7)-24)/4)+8*(row/7);
            new G4LogicalBorderSurface("Grid-AirSurface", short_grid_thick[grid_index], air_inner_x2[iBar], OpGridAirSurface);
            new G4LogicalBorderSurface("Air-GridSurface", air_inner_x2[iBar], short_grid_thick[grid_index], OpAirGridSurface);
        }else{ // short thin grid parts
            grid_index = (44-16*row)*(1-(row/4))+(14+16*(row-4))*(row/4)-4*(iBar-8*row)*(1-((iBar-8*row)/4))-11*((iBar-8*row)/4)+4*(iBar-8*row-4)*((iBar-8*row)/4);
            new G4LogicalBorderSurface("Grid-AirSurface", short_grid_thin[grid_index], air_inner_x2[iBar], OpGridAirSurface);
            new G4LogicalBorderSurface("Air-GridSurface", air_inner_x2[iBar], short_grid_thin[grid_index], OpAirGridSurface);
        }
        if(col==0 || col==4){ // long thick grid parts
            new G4LogicalBorderSurface("Grid-AirSurface", long_grid_thick[col/4], air_inner_y1[iBar], OpGridAirSurface);
            new G4LogicalBorderSurface("Air-GridSurface", air_inner_y1[iBar], long_grid_thick[col/4], OpAirGridSurface);
        }else{ // long thin grid parts
            grid_index = (7-2*col)*(1-col/5)+2*(col-5)*(col/5);
            new G4LogicalBorderSurface("Grid-AirSurface", long_grid_thin[grid_index], air_inner_y1[iBar], OpGridAirSurface);
            new G4LogicalBorderSurface("Air-GridSurface", air_inner_y1[iBar], long_grid_thin[grid_index], OpAirGridSurface);
        }
        if(col==3 || col==7){ // long thick grid parts
            new G4LogicalBorderSurface("Grid-AirSurface", long_grid_thick[(col+1)/4], air_inner_y2[iBar], OpGridAirSurface);
            new G4LogicalBorderSurface("Air-GridSurface", air_inner_y2[iBar], long_grid_thick[(col+1)/4], OpAirGridSurface);
        }else{ // long thin grid parts
            grid_index = (5-2*col)*(1-col/4)+2*(col-4)*(col/4);
            new G4LogicalBorderSurface("Grid-AirSurface", long_grid_thin[grid_index], air_inner_y2[iBar], OpGridAirSurface);
            new G4LogicalBorderSurface("Air-GridSurface", air_inner_y2[iBar], long_grid_thin[grid_index], OpAirGridSurface);
        }
    }
*/
    // ########################################################################################################################################

/*
    // ################# Actual simulations #################

        // Defining arrays
    G4VPhysicalVolume* scintillator_bars[nScintBars];
    G4VPhysicalVolume* sipm_channels[nScintBars];
    G4VPhysicalVolume* vikuiti_x1[nScintBars];
    G4VPhysicalVolume* vikuiti_x2[nScintBars];
    G4VPhysicalVolume* vikuiti_y1[nScintBars];
    G4VPhysicalVolume* vikuiti_y2[nScintBars];
    G4VPhysicalVolume* air_inner_x1[nScintBars];
    G4VPhysicalVolume* air_inner_x2[nScintBars];
    G4VPhysicalVolume* air_inner_y1[nScintBars];
    G4VPhysicalVolume* air_inner_y2[nScintBars];
    G4VPhysicalVolume* air_outer_x1[nScintBars];
    G4VPhysicalVolume* air_outer_x2[nScintBars];
    G4VPhysicalVolume* air_outer_y1[nScintBars];
    G4VPhysicalVolume* air_outer_y2[nScintBars];
    G4VPhysicalVolume* long_grid_thick[3];
    G4VPhysicalVolume* long_grid_thin[6];
    G4VPhysicalVolume* short_grid_thick[24];
    G4VPhysicalVolume* short_grid_thin[48];

        // Scintillator bars, wrapping and sipm
    const G4double x_start = -(BarPitch * 3 + (BarPitch + 0.3) / 2);
    const G4double y_start = -(BarPitch * 3 + (BarPitch + 0.3) / 2);
    for (int iBar = 0; iBar < nScintBars; iBar++) {
        G4double bar_x = x_start + (iBar % 8) * BarPitch;
        G4double bar_y = y_start + (iBar / 8) * BarPitch;
        if (iBar%8 > 3){bar_x += 0.3;}  // the pitch in the middle cross is 6.5mm instead of 6.2mm due to the space of 0.1mm betwenn arrays (plus we have the 200um twice)
        if (iBar/8 > 3){bar_y += 0.3;}
        scintillator_bars[iBar] = new G4PVPlacement(NULL, G4ThreeVector(bar_x, bar_y, 0), ScintLog_, "Scintillator", DetectorLog_, false, iBar);
            // Vikuiti wrapping
        vikuiti_x1[iBar] = new G4PVPlacement(NULL, G4ThreeVector(bar_x, bar_y-(Scint_thick/2.+0.05+Vikuiti_thickness/2.), Grid_height/2.*mm), vikuit_x_log, "Vikuiti", DetectorLog_, false, 2*iBar);
        vikuiti_x2[iBar] = new G4PVPlacement(NULL, G4ThreeVector(bar_x, bar_y+(Scint_thick/2.+0.05+Vikuiti_thickness/2.), Grid_height/2.*mm), vikuit_x_log, "Vikuiti", DetectorLog_, false, 2*iBar+1);
        vikuiti_y1[iBar] = new G4PVPlacement(NULL, G4ThreeVector(bar_x-(Scint_thick/2.+0.05+Vikuiti_thickness/2.), bar_y, Grid_height/2.*mm), vikuit_y_log, "Vikuiti", DetectorLog_, false, 2*iBar);
        vikuiti_y2[iBar] = new G4PVPlacement(NULL, G4ThreeVector(bar_x+(Scint_thick/2.+0.05+Vikuiti_thickness/2.), bar_y, Grid_height/2.*mm), vikuit_y_log, "Vikuiti", DetectorLog_, false, 2*iBar+1);
            // Air strips
        air_inner_x1[iBar] = new G4PVPlacement(NULL, G4ThreeVector(bar_x, bar_y-(Scint_thick/2.+0.05/2.), 0), air_inner_x_log, "Air", DetectorLog_, false, 2*iBar);
        air_inner_x2[iBar] = new G4PVPlacement(NULL, G4ThreeVector(bar_x, bar_y+(Scint_thick/2.+0.05/2.), 0), air_inner_x_log, "Air", DetectorLog_, false, 2*iBar+1);
        air_inner_y1[iBar] = new G4PVPlacement(NULL, G4ThreeVector(bar_x-(Scint_thick/2.+0.05/2.), bar_y, 0), air_inner_y_log, "Air", DetectorLog_, false, 2*iBar);
        air_inner_y2[iBar] = new G4PVPlacement(NULL, G4ThreeVector(bar_x+(Scint_thick/2.+0.05/2.), bar_y, 0), air_inner_y_log, "Air", DetectorLog_, false, 2*iBar+1);
        air_outer_x1[iBar] = new G4PVPlacement(NULL, G4ThreeVector(bar_x, bar_y-(Scint_thick/2.+0.05+Vikuiti_thickness+0.02/2.), 0), air_outer_x_log, "Air", DetectorLog_, false, 2*iBar);
        air_outer_x2[iBar] = new G4PVPlacement(NULL, G4ThreeVector(bar_x, bar_y+(Scint_thick/2.+0.05+Vikuiti_thickness+0.02/2.), 0), air_outer_x_log, "Air", DetectorLog_, false, 2*iBar+1);
        air_outer_y1[iBar] = new G4PVPlacement(NULL, G4ThreeVector(bar_x-(Scint_thick/2.+0.05+Vikuiti_thickness+0.02/2.), bar_y, 0), air_outer_y_log, "Air", DetectorLog_, false, 2*iBar);
        air_outer_y2[iBar] = new G4PVPlacement(NULL, G4ThreeVector(bar_x+(Scint_thick/2.+0.05+Vikuiti_thickness+0.02/2.), bar_y, 0), air_outer_y_log, "Air", DetectorLog_, false, 2*iBar+1);
            // SiPM channels
        sipm_channels[iBar] = new G4PVPlacement(NULL, G4ThreeVector(bar_x, bar_y, -Pad_height-EpoxyResin_height-Scint_length/2.-SiPMsilicon_height/2.*mm), SiPMLog_, "SiPM", DetectorLog_, false, iBar);
    }

        // Vikuiti top layer
    G4VPhysicalVolume* vikuiti_top = new G4PVPlacement(NULL, G4ThreeVector(0, 0, 0.05+(Scint_length+Vikuiti_thickness)/2.*mm), vikuit_z_log, "Vikuiti", DetectorLog_, false, 0);
    G4VPhysicalVolume* air_top = new G4PVPlacement(NULL, G4ThreeVector(0, 0, (Scint_length+0.05)/2.*mm), air_z_log, "Air", DetectorLog_, false, 0);
    G4VPhysicalVolume* air_top2 = new G4PVPlacement(NULL, G4ThreeVector(0, 0, (Scint_length+0.05)/2.+Vikuiti_thickness+0.05*mm), air_z_log, "Air", DetectorLog_, false, 1);


        // Alignement grid
    for (int ii = 0; ii < 3; ii++) {
        long_grid_thick[ii] = new G4PVPlacement(NULL, G4ThreeVector((ii-1)*25.1*mm, 0, -Scint_length/2.+Grid_height/2.*mm), l_grid_long_thick, "Grid", DetectorLog_, false, ii);   
        long_grid_thin[2*ii] = new G4PVPlacement(NULL, G4ThreeVector((6.35+ii*BarPitch)*mm, 0, -Scint_length/2.+Grid_height/2.*mm), l_grid_long_thin, "Grid", DetectorLog_, false, 2*ii);   
        long_grid_thin[2*ii+1] = new G4PVPlacement(NULL, G4ThreeVector(-(6.35+ii*BarPitch)*mm, 0, -Scint_length/2.+Grid_height/2.*mm), l_grid_long_thin, "Grid", DetectorLog_, false, 2*ii+1);   
        for (int jj = 0; jj < 4; jj++) {
            short_grid_thin[16*ii+4*jj] = new G4PVPlacement(NULL, G4ThreeVector(-(3.25+jj*BarPitch)*mm, -(6.35+ii*BarPitch)*mm, -Scint_length/2.+Grid_height/2.*mm), l_grid_short_thin, "Grid", DetectorLog_, false, 16*ii+4*jj);
            short_grid_thin[16*ii+4*jj+1] = new G4PVPlacement(NULL, G4ThreeVector((3.25+jj*BarPitch)*mm, -(6.35+ii*BarPitch)*mm, -Scint_length/2.+Grid_height/2.*mm), l_grid_short_thin, "Grid", DetectorLog_, false, 16*ii+4*jj+1);
            short_grid_thin[16*ii+4*jj+2] = new G4PVPlacement(NULL, G4ThreeVector(-(3.25+jj*BarPitch)*mm, (6.35+ii*BarPitch)*mm, -Scint_length/2.+Grid_height/2.*mm), l_grid_short_thin, "Grid", DetectorLog_, false, 16*ii+4*jj+2);
            short_grid_thin[16*ii+4*jj+3] = new G4PVPlacement(NULL, G4ThreeVector((3.25+jj*BarPitch)*mm, (6.35+ii*BarPitch)*mm, -Scint_length/2.+Grid_height/2.*mm), l_grid_short_thin, "Grid", DetectorLog_, false, 16*ii+4*jj+3);
            short_grid_thick[8*ii+2*jj] = new G4PVPlacement(NULL, G4ThreeVector(-(3.25+jj*BarPitch)*mm, (ii-1)*25.1*mm, -Scint_length/2.+Grid_height/2.*mm), l_grid_short_thick, "Grid", DetectorLog_, false, 8*ii+2*jj);
            short_grid_thick[8*ii+2*jj+1] = new G4PVPlacement(NULL, G4ThreeVector((3.25+jj*BarPitch)*mm, (ii-1)*25.1*mm, -Scint_length/2.+Grid_height/2.*mm), l_grid_short_thick, "Grid", DetectorLog_, false, 8*ii+2*jj+1);
        }
    }

        // Alignement cross
    new G4PVPlacement(NULL, G4ThreeVector(0, 0, Scint_length/2.-Grid_height/2.*mm), l_cross_long, "Cross", DetectorLog_, false, 0);
    new G4PVPlacement(NULL, G4ThreeVector(-45/2/2.*mm, 0, Scint_length/2.-Grid_height/2.*mm), l_cross_short, "Cross", DetectorLog_, false, 0);
    new G4PVPlacement(NULL, G4ThreeVector(45/2/2.*mm, 0, Scint_length/2.-Grid_height/2.*mm), l_cross_short, "Cross", DetectorLog_, false, 1);

        // Optical pad
    G4VPhysicalVolume* optical_pad = new G4PVPlacement(NULL, G4ThreeVector(0,0, -Scint_length/2.-Pad_height/2.*mm), opt_pad_log, "OpticalPad", DetectorLog_, false, 0);
    
        // Epoxy resin SiPM
    G4VPhysicalVolume* epoxy_resin_0 = new G4PVPlacement(NULL, G4ThreeVector(25.1/2.*mm, 25.1/2*mm, -Scint_length/2.-Pad_height-EpoxyResin_height/2.*mm), EpoxyResin_log, "EpoxyResin", DetectorLog_, false, 0);
    G4VPhysicalVolume* epoxy_resin_1 = new G4PVPlacement(NULL, G4ThreeVector(25.1/2.*mm, -25.1/2*mm, -Scint_length/2.-Pad_height-EpoxyResin_height/2.*mm), EpoxyResin_log, "EpoxyResin", DetectorLog_, false, 1);
    G4VPhysicalVolume* epoxy_resin_2 = new G4PVPlacement(NULL, G4ThreeVector(-25.1/2.*mm, 25.1/2*mm, -Scint_length/2.-Pad_height-EpoxyResin_height/2.*mm), EpoxyResin_log, "EpoxyResin", DetectorLog_, false, 2);
    G4VPhysicalVolume* epoxy_resin_3 = new G4PVPlacement(NULL, G4ThreeVector(-25.1/2.*mm, -25.1/2*mm, -Scint_length/2.-Pad_height-EpoxyResin_height/2.*mm), EpoxyResin_log, "EpoxyResin", DetectorLog_, false, 3);


        // Defining interfaces
    // top vikuiti layer
    new G4LogicalBorderSurface("Vikuiti-AirSurface", vikuiti_top, air_top, OpWrapAirSurface);
    new G4LogicalBorderSurface("Air-VikuitiSurface", air_top, vikuiti_top, OpAirWrapSurface);
    new G4LogicalBorderSurface("Vikuiti-AirSurface", vikuiti_top, air_top2, OpWrapAirSurface);
    new G4LogicalBorderSurface("Air-VikuitiSurface", air_top2, vikuiti_top, OpAirWrapSurface);
    // optical pad - resin
    new G4LogicalBorderSurface("ResinPadSurface", epoxy_resin_0, optical_pad, OpResinPadSurface);
    new G4LogicalBorderSurface("PadResinSurface", optical_pad, epoxy_resin_0, OpPadResinSurface);
    new G4LogicalBorderSurface("ResinPadSurface", epoxy_resin_1, optical_pad, OpResinPadSurface);
    new G4LogicalBorderSurface("PadResinSurface", optical_pad, epoxy_resin_1, OpPadResinSurface);
    new G4LogicalBorderSurface("ResinPadSurface", epoxy_resin_2, optical_pad, OpResinPadSurface);
    new G4LogicalBorderSurface("PadResinSurface", optical_pad, epoxy_resin_2, OpPadResinSurface);
    new G4LogicalBorderSurface("ResinPadSurface", epoxy_resin_3, optical_pad, OpResinPadSurface);
    new G4LogicalBorderSurface("PadResinSurface", optical_pad, epoxy_resin_3, OpPadResinSurface);
    for (int iBar = 0; iBar < nScintBars; iBar++) {
        G4int col = iBar%8;
        G4int row = iBar/8; 
        // scintillator 4 sides
        new G4LogicalBorderSurface("Sci-AirSurface", scintillator_bars[iBar], air_inner_x1[iBar], OpSciAirSurface);
        new G4LogicalBorderSurface("Air-SciSurface", air_inner_x1[iBar], scintillator_bars[iBar], OpAirSciSurface);
        new G4LogicalBorderSurface("Sci-AirSurface", scintillator_bars[iBar], air_inner_x2[iBar], OpSciAirSurface);
        new G4LogicalBorderSurface("Air-SciSurface", air_inner_x2[iBar], scintillator_bars[iBar], OpAirSciSurface);
        new G4LogicalBorderSurface("Sci-AirSurface", scintillator_bars[iBar], air_inner_y1[iBar], OpSciAirSurface);
        new G4LogicalBorderSurface("Air-SciSurface", air_inner_y1[iBar], scintillator_bars[iBar], OpAirSciSurface);
        new G4LogicalBorderSurface("Sci-AirSurface", scintillator_bars[iBar], air_inner_y2[iBar], OpSciAirSurface);
        new G4LogicalBorderSurface("Air-SciSurface", air_inner_y2[iBar], scintillator_bars[iBar], OpAirSciSurface);
        // scintillator top
        new G4LogicalBorderSurface("Sci-AirSurface", scintillator_bars[iBar], air_top, OpSciAirSurface);
        new G4LogicalBorderSurface("Air-SciSurface", air_top, scintillator_bars[iBar], OpAirSciSurface);
        // scintillator bottom (with optical pad)
        new G4LogicalBorderSurface("Sci-PadSurface", scintillator_bars[iBar], optical_pad, OpSciPadSurface);
        new G4LogicalBorderSurface("Pad-SciSurface", optical_pad, scintillator_bars[iBar], OpPadSciSurface);
        // wrapping
        new G4LogicalBorderSurface("Vikuiti-AirSurface", vikuiti_x1[iBar], air_inner_x1[iBar], OpWrapAirSurface);
        new G4LogicalBorderSurface("Vikuiti-AirSurface", vikuiti_x1[iBar], air_outer_x1[iBar], OpWrapAirSurface);
        new G4LogicalBorderSurface("Air-VikuitiSurface", air_inner_x1[iBar], vikuiti_x1[iBar], OpAirWrapSurface);
        new G4LogicalBorderSurface("Air-VikuitiSurface", air_outer_x1[iBar], vikuiti_x1[iBar], OpAirWrapSurface);
        new G4LogicalBorderSurface("Vikuiti-AirSurface", vikuiti_x2[iBar], air_inner_x2[iBar], OpWrapAirSurface);
        new G4LogicalBorderSurface("Vikuiti-AirSurface", vikuiti_x2[iBar], air_outer_x2[iBar], OpWrapAirSurface);
        new G4LogicalBorderSurface("Air-VikuitiSurface", air_inner_x2[iBar], vikuiti_x2[iBar], OpAirWrapSurface);
        new G4LogicalBorderSurface("Air-VikuitiSurface", air_outer_x2[iBar], vikuiti_x2[iBar], OpAirWrapSurface);
        new G4LogicalBorderSurface("Vikuiti-AirSurface", vikuiti_y1[iBar], air_inner_y1[iBar], OpWrapAirSurface);
        new G4LogicalBorderSurface("Vikuiti-AirSurface", vikuiti_y1[iBar], air_outer_y1[iBar], OpWrapAirSurface);
        new G4LogicalBorderSurface("Air-VikuitiSurface", air_inner_y1[iBar], vikuiti_y1[iBar], OpAirWrapSurface);
        new G4LogicalBorderSurface("Air-VikuitiSurface", air_outer_y1[iBar], vikuiti_y1[iBar], OpAirWrapSurface);
        new G4LogicalBorderSurface("Vikuiti-AirSurface", vikuiti_y2[iBar], air_inner_y2[iBar], OpWrapAirSurface);
        new G4LogicalBorderSurface("Vikuiti-AirSurface", vikuiti_y2[iBar], air_outer_y2[iBar], OpWrapAirSurface);
        new G4LogicalBorderSurface("Air-VikuitiSurface", air_inner_y2[iBar], vikuiti_y2[iBar], OpAirWrapSurface);
        new G4LogicalBorderSurface("Air-VikuitiSurface", air_outer_y2[iBar], vikuiti_y2[iBar], OpAirWrapSurface);
        // sipm PDE (b/w SiPM and resin)
        if (row<4){
          if (col<4){
            new G4LogicalBorderSurface("ResinSiPMSurface", epoxy_resin_3, sipm_channels[iBar], OpResinSiPMSurface);
            new G4LogicalBorderSurface("ResinSiPMSurface", sipm_channels[iBar], epoxy_resin_3, OpSiPMResinSurface);
          }else{
            new G4LogicalBorderSurface("ResinSiPMSurface", epoxy_resin_1, sipm_channels[iBar], OpResinSiPMSurface);
            new G4LogicalBorderSurface("ResinSiPMSurface", sipm_channels[iBar], epoxy_resin_1, OpSiPMResinSurface);
          }
        }else{
          if (col<4){
            new G4LogicalBorderSurface("ResinSiPMSurface", epoxy_resin_2, sipm_channels[iBar], OpResinSiPMSurface);
            new G4LogicalBorderSurface("ResinSiPMSurface", sipm_channels[iBar], epoxy_resin_2, OpSiPMResinSurface);
          }else{
            new G4LogicalBorderSurface("ResinSiPMSurface", epoxy_resin_0, sipm_channels[iBar], OpResinSiPMSurface);
            new G4LogicalBorderSurface("ResinSiPMSurface", sipm_channels[iBar], epoxy_resin_0, OpSiPMResinSurface);
          }
        }
*/        /*if (row<4){
          if (col<4){
            new G4LogicalBorderSurface("ResinSiPMSurface", scintillator_bars[iBar], sipm_channels[iBar], OpResinSiPMSurface);
            new G4LogicalBorderSurface("ResinSiPMSurface", sipm_channels[iBar], scintillator_bars[iBar], OpSiPMResinSurface);
          }else{
            new G4LogicalBorderSurface("ResinSiPMSurface", scintillator_bars[iBar], sipm_channels[iBar], OpResinSiPMSurface);
            new G4LogicalBorderSurface("ResinSiPMSurface", sipm_channels[iBar], scintillator_bars[iBar], OpSiPMResinSurface);
          }
        }else{
          if (col<4){
            new G4LogicalBorderSurface("ResinSiPMSurface", scintillator_bars[iBar], sipm_channels[iBar], OpResinSiPMSurface);
            new G4LogicalBorderSurface("ResinSiPMSurface", sipm_channels[iBar], scintillator_bars[iBar], OpSiPMResinSurface);
          }else{
            new G4LogicalBorderSurface("ResinSiPMSurface", scintillator_bars[iBar], sipm_channels[iBar], OpResinSiPMSurface);
            new G4LogicalBorderSurface("ResinSiPMSurface", sipm_channels[iBar], scintillator_bars[iBar], OpSiPMResinSurface);
          }
        }*/
/*        // grid
        G4int grid_index;
        if(row==0 || row==4){ // short thick grid parts
            grid_index = (6-2*(iBar-32*(row/4)))*(1-(iBar-32*(row/4))/4)+(2*(iBar-32*(row/4)-4)+1)*((iBar-32*(row/4))/4)+8*(row/4);
            new G4LogicalBorderSurface("Grid-AirSurface", short_grid_thick[grid_index], air_inner_x1[iBar], OpGridAirSurface);
            new G4LogicalBorderSurface("Air-GridSurface", air_inner_x1[iBar], short_grid_thick[grid_index], OpAirGridSurface);
        }else{ // short thin grid parts
            grid_index = (44-16*(row-1))*(1-(row/4))+(14+16*(row-5))*(row/4)-4*(iBar-8*row)*(1-((iBar-8*row)/4))-11*((iBar-8*row)/4)+4*(iBar-8*row-4)*((iBar-8*row)/4);
            new G4LogicalBorderSurface("Grid-AirSurface", short_grid_thin[grid_index], air_inner_x1[iBar], OpGridAirSurface);
            new G4LogicalBorderSurface("Air-GridSurface", air_inner_x1[iBar], short_grid_thin[grid_index], OpAirGridSurface);
        }
        if(row==3 || row==7){ // short thick grid parts
            grid_index = (14-2*(iBar-32*(row/7)-24))*(1-((iBar-32*(row/7)-24)/4))+(2*(iBar-32*(row/7)-28)+9)*((iBar-32*(row/7)-24)/4)+8*(row/7);
            new G4LogicalBorderSurface("Grid-AirSurface", short_grid_thick[grid_index], air_inner_x2[iBar], OpGridAirSurface);
            new G4LogicalBorderSurface("Air-GridSurface", air_inner_x2[iBar], short_grid_thick[grid_index], OpAirGridSurface);
        }else{ // short thin grid parts
            grid_index = (44-16*row)*(1-(row/4))+(14+16*(row-4))*(row/4)-4*(iBar-8*row)*(1-((iBar-8*row)/4))-11*((iBar-8*row)/4)+4*(iBar-8*row-4)*((iBar-8*row)/4);
            new G4LogicalBorderSurface("Grid-AirSurface", short_grid_thin[grid_index], air_inner_x2[iBar], OpGridAirSurface);
            new G4LogicalBorderSurface("Air-GridSurface", air_inner_x2[iBar], short_grid_thin[grid_index], OpAirGridSurface);
        }
        if(col==0 || col==4){ // long thick grid parts
            new G4LogicalBorderSurface("Grid-AirSurface", long_grid_thick[col/4], air_inner_y1[iBar], OpGridAirSurface);
            new G4LogicalBorderSurface("Air-GridSurface", air_inner_y1[iBar], long_grid_thick[col/4], OpAirGridSurface);
        }else{ // long thin grid parts
            grid_index = (7-2*col)*(1-col/5)+2*(col-5)*(col/5);
            new G4LogicalBorderSurface("Grid-AirSurface", long_grid_thin[grid_index], air_inner_y1[iBar], OpGridAirSurface);
            new G4LogicalBorderSurface("Air-GridSurface", air_inner_y1[iBar], long_grid_thin[grid_index], OpAirGridSurface);
        }
        if(col==3 || col==7){ // long thick grid parts
            new G4LogicalBorderSurface("Grid-AirSurface", long_grid_thick[(col+1)/4], air_inner_y2[iBar], OpGridAirSurface);
            new G4LogicalBorderSurface("Air-GridSurface", air_inner_y2[iBar], long_grid_thick[(col+1)/4], OpAirGridSurface);
        }else{ // long thin grid parts
            grid_index = (5-2*col)*(1-col/4)+2*(col-4)*(col/4);
            new G4LogicalBorderSurface("Grid-AirSurface", long_grid_thin[grid_index], air_inner_y2[iBar], OpGridAirSurface);
            new G4LogicalBorderSurface("Air-GridSurface", air_inner_y2[iBar], long_grid_thin[grid_index], OpAirGridSurface);
        }
    }



    // define logical surfaces b/w SiPM resin/optical pad/scintillator ??

*/


    

/*
    // ################# perfect detector (perfect reflector, no grid, no crosstalk, SiPM directly in contact with bars) #################
        // Defining arrays
    G4VPhysicalVolume* scintillator_bars[nScintBars];
    G4VPhysicalVolume* PerfRefl_x1[nScintBars];
    G4VPhysicalVolume* PerfRefl_x2[nScintBars];
    G4VPhysicalVolume* PerfRefl_y1[nScintBars];
    G4VPhysicalVolume* PerfRefl_y2[nScintBars];
    G4VPhysicalVolume* air_inner_x1[nScintBars];
    G4VPhysicalVolume* air_inner_x2[nScintBars];
    G4VPhysicalVolume* air_inner_y1[nScintBars];
    G4VPhysicalVolume* air_inner_y2[nScintBars];
    G4VPhysicalVolume* air_outer_x1[nScintBars];
    G4VPhysicalVolume* air_outer_x2[nScintBars];
    G4VPhysicalVolume* air_outer_y1[nScintBars];
    G4VPhysicalVolume* air_outer_y2[nScintBars];


        // Scintillator bars and wrapping
    const G4double x_start = -(BarPitch * 3 + (BarPitch + 0.3) / 2);
    const G4double y_start = -(BarPitch * 3 + (BarPitch + 0.3) / 2);
    for (int iBar = 0; iBar < nScintBars; iBar++) {
        G4double bar_x = x_start + (iBar % 8) * BarPitch;
        G4double bar_y = y_start + (iBar / 8) * BarPitch;
        if (iBar%8 > 3){bar_x += 0.3;}  // the pitch in the middle cross is 6.5mm instead of 6.2mm due to the space of 0.1mm betwenn arrays (plus we have the 200um twice)
        if (iBar/8 > 3){bar_y += 0.3;}
        scintillator_bars[iBar] = new G4PVPlacement(NULL, G4ThreeVector(bar_x, bar_y, 0), ScintLog_, "Scintillator", DetectorLog_, false, iBar);
            // Perfect reflector wrapping
        PerfRefl_x1[iBar] = new G4PVPlacement(NULL, G4ThreeVector(bar_x, bar_y-(Scint_thick/2.+0.05+Vikuiti_thickness/2.), 0), PerfRefl_x_log, "PerfRefl", DetectorLog_, false, 2*iBar);
        PerfRefl_x2[iBar] = new G4PVPlacement(NULL, G4ThreeVector(bar_x, bar_y+(Scint_thick/2.+0.05+Vikuiti_thickness/2.), 0), PerfRefl_x_log, "PerfRefl", DetectorLog_, false, 2*iBar+1);
        PerfRefl_y1[iBar] = new G4PVPlacement(NULL, G4ThreeVector(bar_x-(Scint_thick/2.+0.05+Vikuiti_thickness/2.), bar_y, 0), PerfRefl_y_log, "PerfRefl", DetectorLog_, false, 2*iBar);
        PerfRefl_y2[iBar] = new G4PVPlacement(NULL, G4ThreeVector(bar_x+(Scint_thick/2.+0.05+Vikuiti_thickness/2.), bar_y, 0), PerfRefl_y_log, "PerfRefl", DetectorLog_, false, 2*iBar+1);
            // Air strips
        air_inner_x1[iBar] = new G4PVPlacement(NULL, G4ThreeVector(bar_x, bar_y-(Scint_thick/2.+0.05/2.), 0), air_inner_x_log, "Air", DetectorLog_, false, 2*iBar);
        air_inner_x2[iBar] = new G4PVPlacement(NULL, G4ThreeVector(bar_x, bar_y+(Scint_thick/2.+0.05/2.), 0), air_inner_x_log, "Air", DetectorLog_, false, 2*iBar+1);
        air_inner_y1[iBar] = new G4PVPlacement(NULL, G4ThreeVector(bar_x-(Scint_thick/2.+0.05/2.), bar_y, 0), air_inner_y_log, "Air", DetectorLog_, false, 2*iBar);
        air_inner_y2[iBar] = new G4PVPlacement(NULL, G4ThreeVector(bar_x+(Scint_thick/2.+0.05/2.), bar_y, 0), air_inner_y_log, "Air", DetectorLog_, false, 2*iBar+1);
        air_outer_x1[iBar] = new G4PVPlacement(NULL, G4ThreeVector(bar_x, bar_y-(Scint_thick/2.+0.05+Vikuiti_thickness+0.02/2.), 0), air_outer_x_log, "Air", DetectorLog_, false, 2*iBar);
        air_outer_x2[iBar] = new G4PVPlacement(NULL, G4ThreeVector(bar_x, bar_y+(Scint_thick/2.+0.05+Vikuiti_thickness+0.02/2.), 0), air_outer_x_log, "Air", DetectorLog_, false, 2*iBar+1);
        air_outer_y1[iBar] = new G4PVPlacement(NULL, G4ThreeVector(bar_x-(Scint_thick/2.+0.05+Vikuiti_thickness+0.02/2.), bar_y, 0), air_outer_y_log, "Air", DetectorLog_, false, 2*iBar);
        air_outer_y2[iBar] = new G4PVPlacement(NULL, G4ThreeVector(bar_x+(Scint_thick/2.+0.05+Vikuiti_thickness+0.02/2.), bar_y, 0), air_outer_y_log, "Air", DetectorLog_, false, 2*iBar+1);
            // SiPM channels
        new G4PVPlacement(NULL, G4ThreeVector(bar_x, bar_y, -Scint_length/2.-SiPMsilicon_height/2.*mm), SiPMLog_, "SiPM", DetectorLog_, false, iBar);
    }

        // PerfRefl top layer
    G4VPhysicalVolume* PerfRefl_top = new G4PVPlacement(NULL, G4ThreeVector(0, 0, 0.05+(Scint_length+Vikuiti_thickness)/2.*mm), PerfRefl_z_log, "PerfRefl", DetectorLog_, false, 0);
    G4VPhysicalVolume* air_top = new G4PVPlacement(NULL, G4ThreeVector(0, 0, (Scint_length+0.05)/2.*mm), air_z_log, "Air", DetectorLog_, false, 0);
    G4VPhysicalVolume* air_top2 = new G4PVPlacement(NULL, G4ThreeVector(0, 0, (Scint_length+0.05)/2.+Vikuiti_thickness+0.05*mm), air_z_log, "Air", DetectorLog_, false, 1);

        // Defining interfaces
    // top PerfRefl layer
    new G4LogicalBorderSurface("PerfRefl-AirSurface", PerfRefl_top, air_top, OpPerfReflAirSurface);
    new G4LogicalBorderSurface("Air-PerfReflSurface", air_top, PerfRefl_top, OpAirPerfReflSurface);
    new G4LogicalBorderSurface("PerfRefl-AirSurface", PerfRefl_top, air_top2, OpPerfReflAirSurface);
    new G4LogicalBorderSurface("Air-PerfReflSurface", air_top2, PerfRefl_top, OpAirPerfReflSurface);
    for (int iBar = 0; iBar < nScintBars; iBar++) {
        // scintillator 4 sides
        new G4LogicalBorderSurface("Sci-AirSurface", scintillator_bars[iBar], air_inner_x1[iBar], OpSciAirSurface);
        new G4LogicalBorderSurface("Air-SciSurface", air_inner_x1[iBar], scintillator_bars[iBar], OpAirSciSurface);
        new G4LogicalBorderSurface("Sci-AirSurface", scintillator_bars[iBar], air_inner_x2[iBar], OpSciAirSurface);
        new G4LogicalBorderSurface("Air-SciSurface", air_inner_x2[iBar], scintillator_bars[iBar], OpAirSciSurface);
        new G4LogicalBorderSurface("Sci-AirSurface", scintillator_bars[iBar], air_inner_y1[iBar], OpSciAirSurface);
        new G4LogicalBorderSurface("Air-SciSurface", air_inner_y1[iBar], scintillator_bars[iBar], OpAirSciSurface);
        new G4LogicalBorderSurface("Sci-AirSurface", scintillator_bars[iBar], air_inner_y2[iBar], OpSciAirSurface);
        new G4LogicalBorderSurface("Air-SciSurface", air_inner_y2[iBar], scintillator_bars[iBar], OpAirSciSurface);
        // scintillator top
        new G4LogicalBorderSurface("Sci-AirSurface", scintillator_bars[iBar], air_top, OpSciAirSurface);
        new G4LogicalBorderSurface("Air-SciSurface", air_top, scintillator_bars[iBar], OpAirSciSurface);
        // wrapping
        new G4LogicalBorderSurface("PerfRefl-AirSurface", PerfRefl_x1[iBar], air_inner_x1[iBar], OpPerfReflAirSurface);
        new G4LogicalBorderSurface("PerfRefl-AirSurface", PerfRefl_x1[iBar], air_outer_x1[iBar], OpPerfReflAirSurface);
        new G4LogicalBorderSurface("Air-PerfReflSurface", air_inner_x1[iBar], PerfRefl_x1[iBar], OpAirPerfReflSurface);
        new G4LogicalBorderSurface("Air-PerfReflSurface", air_outer_x1[iBar], PerfRefl_x1[iBar], OpAirPerfReflSurface);
        new G4LogicalBorderSurface("PerfRefl-AirSurface", PerfRefl_x2[iBar], air_inner_x2[iBar], OpPerfReflAirSurface);
        new G4LogicalBorderSurface("PerfRefl-AirSurface", PerfRefl_x2[iBar], air_outer_x2[iBar], OpPerfReflAirSurface);
        new G4LogicalBorderSurface("Air-PerfReflSurface", air_inner_x2[iBar], PerfRefl_x2[iBar], OpAirPerfReflSurface);
        new G4LogicalBorderSurface("Air-PerfReflSurface", air_outer_x2[iBar], PerfRefl_x2[iBar], OpAirPerfReflSurface);
        new G4LogicalBorderSurface("PerfRefl-AirSurface", PerfRefl_y1[iBar], air_inner_y1[iBar], OpPerfReflAirSurface);
        new G4LogicalBorderSurface("PerfRefl-AirSurface", PerfRefl_y1[iBar], air_outer_y1[iBar], OpPerfReflAirSurface);
        new G4LogicalBorderSurface("Air-PerfReflSurface", air_inner_y1[iBar], PerfRefl_y1[iBar], OpAirPerfReflSurface);
        new G4LogicalBorderSurface("Air-PerfReflSurface", air_outer_y1[iBar], PerfRefl_y1[iBar], OpAirPerfReflSurface);
        new G4LogicalBorderSurface("PerfRefl-AirSurface", PerfRefl_y2[iBar], air_inner_y2[iBar], OpPerfReflAirSurface);
        new G4LogicalBorderSurface("PerfRefl-AirSurface", PerfRefl_y2[iBar], air_outer_y2[iBar], OpPerfReflAirSurface);
        new G4LogicalBorderSurface("Air-PerfReflSurface", air_inner_y2[iBar], PerfRefl_y2[iBar], OpAirPerfReflSurface);
        new G4LogicalBorderSurface("Air-PerfReflSurface", air_outer_y2[iBar], PerfRefl_y2[iBar], OpAirPerfReflSurface);
    }

    // ########################################################################################################################################
*/

/*
    // ################################ DEBUGGING ################################

      // Testing the scintillation efficiency
    new G4PVPlacement(NULL, G4ThreeVector(0, 0, 0), ScintBlockLog_, "Scintillator", DetectorLog_, false, 0);


      // Testing the SiPM PDE
    G4VPhysicalVolume* theSiPM = new G4PVPlacement(NULL, G4ThreeVector(0,0, -20*mm), SiPMLog_, "SiPM", DetectorLog_, false, 0);
    G4VPhysicalVolume* theAir = new G4PVPlacement(NULL, G4ThreeVector(0, 0, -20+SiPMsilicon_height/2+2/2.*mm), air_test_log, "Air", DetectorLog_, false, 0);

    new G4LogicalBorderSurface("SiPMSurface", theAir, theSiPM, OpResinSiPMSurface);


      // Checking the Perfect Reflector reflectivity/transmittance
    G4VPhysicalVolume* perfectrefl = new G4PVPlacement(NULL, G4ThreeVector(0, 0, 0), PerfRefl_z_log, "PerfRefl", DetectorLog_, false, 0);

    G4VPhysicalVolume* airtop = new G4PVPlacement(NULL, G4ThreeVector(0, 0, 1+Vikuiti_thickness/2.*mm), air_test_log, "Air", DetectorLog_, false, 0);
    G4VPhysicalVolume* airbottom = new G4PVPlacement(NULL, G4ThreeVector(0, 0, -1-Vikuiti_thickness/2.*mm), air_test_log, "Air", DetectorLog_, false, 1);

    new G4LogicalBorderSurface("Air-PerfReflSurface", airtop, perfectrefl, OpAirPerfReflSurface);
    new G4LogicalBorderSurface("Air-PerfReflSurface", airbottom, perfectrefl, OpAirPerfReflSurface);
    new G4LogicalBorderSurface("PerfRefl-AirSurface", perfectrefl, airtop, OpPerfReflAirSurface);
    new G4LogicalBorderSurface("PerfRefl-AirSurface", perfectrefl, airbottom, OpPerfReflAirSurface);

    new G4PVPlacement(NULL, G4ThreeVector(0,0, -20*mm), SiPMLog_, "SiPM", DetectorLog_, false, 0);
    new G4PVPlacement(NULL, G4ThreeVector(0,0, 20*mm), SiPMLog_, "SiPM", DetectorLog_, false, 1);


      // Tuning the Vikuiti reflectivity/transmittance
    G4VPhysicalVolume* vikuiti = new G4PVPlacement(NULL, G4ThreeVector(0, 0, 0), vikuit_z_log, "Vikuiti", DetectorLog_, false, 0);

    G4VPhysicalVolume* airtop = new G4PVPlacement(NULL, G4ThreeVector(0, 0, 1+Vikuiti_thickness/2.*mm), air_test_log, "Air", DetectorLog_, false, 0);
    G4VPhysicalVolume* airbottom = new G4PVPlacement(NULL, G4ThreeVector(0, 0, -1-Vikuiti_thickness/2.*mm), air_test_log, "Air", DetectorLog_, false, 1);

    new G4LogicalBorderSurface("Air-VikuitiSurface", airtop, vikuiti, OpAirWrapSurface);
    new G4LogicalBorderSurface("Air-VikuitiSurface", airbottom, vikuiti, OpAirWrapSurface);
    new G4LogicalBorderSurface("Vikuiti-AirSurface", vikuiti, airtop, OpWrapAirSurface);
    new G4LogicalBorderSurface("Vikuiti-AirSurface", vikuiti, airbottom, OpWrapAirSurface);

    new G4PVPlacement(NULL, G4ThreeVector(0,0, -20*mm), SiPMLog_, "SiPM", DetectorLog_, false, 0);
    new G4PVPlacement(NULL, G4ThreeVector(0,0, 20*mm), SiPMLog_, "SiPM", DetectorLog_, false, 1);


      // Tuning the Grid transmittance
    G4VPhysicalVolume* thegrid = new G4PVPlacement(NULL, G4ThreeVector(0, 0, 0), grid_plate_test, "Grid", DetectorLog_, false, 0);

    G4VPhysicalVolume* airtop = new G4PVPlacement(NULL, G4ThreeVector(0, 0, 1+0.2/2.*mm), air_test_log, "Air", DetectorLog_, false, 0);
    G4VPhysicalVolume* airbottom = new G4PVPlacement(NULL, G4ThreeVector(0, 0, -1-0.2/2.*mm), air_test_log, "Air", DetectorLog_, false, 1);

    new G4LogicalBorderSurface("Air-GridSurface", airtop, thegrid, OpAirGridSurface);
    new G4LogicalBorderSurface("Air-GridSurface", airbottom, thegrid, OpAirGridSurface);
    new G4LogicalBorderSurface("Grid-AirSurface", thegrid, airtop, OpGridAirSurface);
    new G4LogicalBorderSurface("Grid-AirSurface", thegrid, airbottom, OpGridAirSurface);

    new G4PVPlacement(NULL, G4ThreeVector(0,0, -20*mm), SiPMLog_, "SiPM", DetectorLog_, false, 0);
    new G4PVPlacement(NULL, G4ThreeVector(0,0, 20*mm), SiPMLog_, "SiPM", DetectorLog_, false, 1);
*/


    // ################################ LEAP SCINTILLATOR ################################

    /*G4int xy_shift = 0;
    G4int display_split = 0;
    if (LEAPvisualization){  // for visualization purposes
      xy_shift = 200;
      display_split = 20;
    }

    G4VPhysicalVolume* LEAPscintillator_bar = new G4PVPlacement(NULL, G4ThreeVector(xy_shift, xy_shift, 0), LEAPScintLog_, "LEAPScintillator", DetectorLog_, false, 0);

    G4VPhysicalVolume* LEAP_Optpad = new G4PVPlacement(NULL, G4ThreeVector(xy_shift, xy_shift, -display_split-LEAPScint_length/2.-LEAP_Optpad_thick/2.), LEAPopt_pad_log, "Optpad", DetectorLog_, false, 0);
    G4VPhysicalVolume* LEAP_PMT = new G4PVPlacement(NULL, G4ThreeVector(xy_shift, xy_shift, -3*display_split-LEAPScint_length/2.-LEAP_Optpad_thick-LEAP_PMT_thick/2.), LEAPPMTLog_, "PMT", DetectorLog_, false, 0);
 
    G4double Teflon_shift;
    if (LEAPScint_conical){
      Teflon_shift = -2*display_split-LEAPScint_length/2.+LEAPTefloncone_thick/2.-LEAPTeflonsquare_thick-0.05;
    } else{
      Teflon_shift = -2*display_split-LEAPScint_length/2.-LEAPTeflonsquare_thick/2.-0.05;
    }
    G4VPhysicalVolume* LEAP_TeflonAdapter = new G4PVPlacement(NULL, G4ThreeVector(xy_shift, xy_shift, Teflon_shift), TeflonLog_, "TeflonAdapter", DetectorLog_, false, 0);
    
    G4double air_teflon_shift;
    if (LEAPScint_conical){
      air_teflon_shift = -LEAPScint_length/2.-0.05/2.+LEAPTefloncone_thick-LEAPTeflonsquare_thick;
    } else{
      air_teflon_shift = -LEAPScint_length/2.-0.05/2.;
    }
    G4VPhysicalVolume* air_teflon = new G4PVPlacement(NULL, G4ThreeVector(xy_shift, xy_shift, air_teflon_shift), LEAPair_teflonsquare_log, "Air", DetectorLog_, false, 0);
    G4double air_teflon_coneshift;
    G4VPhysicalVolume* air_teflon_cone;
    if (LEAPScint_conical){
      air_teflon_coneshift = -LEAPScint_length/2. + sqrt(3)*(LEAPScint_thick/sqrt(2)-4)/2.;
      air_teflon_cone = new G4PVPlacement(NULL, G4ThreeVector(xy_shift, xy_shift, air_teflon_coneshift), LEAPair_teflonconical_log, "Air", DetectorLog_, false, 0);
    }
    //if (LEAPScint_conical){
      //air_teflon_shift = -LEAPScint_length/2.-0.05/2.+(LEAPTefloncone_thick-LEAPTeflonsquare_thick)/2.;
    //  air_teflon_shift = -LEAPScint_length/2. + sqrt(3)*(LEAPScint_thick/sqrt(2)-4)/2.;
    //  air_teflon = new G4PVPlacement(NULL, G4ThreeVector(xy_shift, xy_shift, air_teflon_shift), LEAPair_teflonconical_log, "Air", DetectorLog_, false, 0);
    //} else{
    //  air_teflon_shift = -LEAPScint_length/2.-0.05/2.;
    //  air_teflon = new G4PVPlacement(NULL, G4ThreeVector(xy_shift, xy_shift, air_teflon_shift), LEAPair_teflonsquare_log, "Air", DetectorLog_, false, 0);
    //}


    G4double vikuiti_shift;
    if (LEAPScint_conical){
      vikuiti_shift = (LEAPTefloncone_thick-LEAPTeflonsquare_thick)/2.;
    } else{
      vikuiti_shift = 0.0*mm;
    }
    G4VPhysicalVolume* LEAPvikuiti_x1 = new G4PVPlacement(NULL, G4ThreeVector(xy_shift, -display_split+xy_shift-(LEAPScint_thick/2.+0.05+Vikuiti_thickness/2.), vikuiti_shift), LEAPvikuit_x_log, "Vikuiti", DetectorLog_, false, 0);
    G4VPhysicalVolume* LEAPvikuiti_x2 = new G4PVPlacement(NULL, G4ThreeVector(xy_shift, display_split+xy_shift+(LEAPScint_thick/2.+0.05+Vikuiti_thickness/2.), vikuiti_shift), LEAPvikuit_x_log, "Vikuiti", DetectorLog_, false, 1);
    G4VPhysicalVolume* LEAPvikuiti_y1 = new G4PVPlacement(NULL, G4ThreeVector(-display_split+xy_shift-(LEAPScint_thick/2.+0.05+Vikuiti_thickness/2.), xy_shift, vikuiti_shift), LEAPvikuit_y_log, "Vikuiti", DetectorLog_, false, 0);
    G4VPhysicalVolume* LEAPvikuiti_y2 = new G4PVPlacement(NULL, G4ThreeVector(display_split+xy_shift+(LEAPScint_thick/2.+0.05+Vikuiti_thickness/2.), xy_shift, vikuiti_shift), LEAPvikuit_y_log, "Vikuiti", DetectorLog_, false, 1);

    G4VPhysicalVolume* LEAPvikuiti_top = new G4PVPlacement(NULL, G4ThreeVector(xy_shift, xy_shift, display_split+0.05+(LEAPScint_length+Vikuiti_thickness)/2.*mm), LEAPvikuit_z_log, "Vikuiti", DetectorLog_, false, 0);
  
    G4double air_shift;
    if (LEAPScint_conical){
      air_shift = (LEAPTefloncone_thick-LEAPTeflonsquare_thick)/2.*mm;
    } else{
      air_shift = 0.0*mm;
    }
    G4VPhysicalVolume* air_x1 = new G4PVPlacement(NULL, G4ThreeVector(xy_shift, xy_shift-(LEAPScint_thick/2.+0.05/2.), air_shift), LEAPair_x_log, "Air", DetectorLog_, false, 0);
    G4VPhysicalVolume* air_x2 = new G4PVPlacement(NULL, G4ThreeVector(xy_shift, xy_shift+(LEAPScint_thick/2.+0.05/2.), air_shift), LEAPair_x_log, "Air", DetectorLog_, false, 1);
    G4VPhysicalVolume* air_y1 = new G4PVPlacement(NULL, G4ThreeVector(xy_shift-(LEAPScint_thick/2.+0.05/2.), xy_shift, air_shift), LEAPair_y_log, "Air", DetectorLog_, false, 0);
    G4VPhysicalVolume* air_y2 = new G4PVPlacement(NULL, G4ThreeVector(xy_shift+(LEAPScint_thick/2.+0.05/2.), xy_shift, air_shift), LEAPair_y_log, "Air", DetectorLog_, false, 1);
    G4VPhysicalVolume* air_top = new G4PVPlacement(NULL, G4ThreeVector(xy_shift, xy_shift, (LEAPScint_length+0.05)/2.*mm), LEAPair_z_log, "Air", DetectorLog_, false, 0);

    // PMT
    new G4LogicalBorderSurface("ResinSiPMSurface", LEAPscintillator_bar, LEAP_PMT, OpResinSiPMSurface);
    new G4LogicalBorderSurface("ResinSiPMSurface", LEAP_PMT, LEAPscintillator_bar, OpSiPMResinSurface);

    // Vikuiti
    new G4LogicalBorderSurface("Vikuiti-AirSurface", LEAPvikuiti_x1, air_x1, OpWrapAirSurface);
    new G4LogicalBorderSurface("Air-VikuitiSurface", air_x1, LEAPvikuiti_x1, OpAirWrapSurface);
    new G4LogicalBorderSurface("Vikuiti-AirSurface", LEAPvikuiti_x2, air_x2, OpWrapAirSurface);
    new G4LogicalBorderSurface("Air-VikuitiSurface", air_x2, LEAPvikuiti_x2, OpAirWrapSurface);
    new G4LogicalBorderSurface("Vikuiti-AirSurface", LEAPvikuiti_y1, air_y1, OpWrapAirSurface);
    new G4LogicalBorderSurface("Air-VikuitiSurface", air_y1, LEAPvikuiti_y1, OpAirWrapSurface);
    new G4LogicalBorderSurface("Vikuiti-AirSurface", LEAPvikuiti_y2, air_y2, OpWrapAirSurface);
    new G4LogicalBorderSurface("Air-VikuitiSurface", air_y2, LEAPvikuiti_y2, OpAirWrapSurface);
    new G4LogicalBorderSurface("Vikuiti-AirSurface", LEAPvikuiti_top, air_top, OpWrapAirSurface);
    new G4LogicalBorderSurface("Air-VikuitiSurface", air_top, LEAPvikuiti_top, OpAirWrapSurface);

    // Scintillator
    new G4LogicalBorderSurface("Sci-AirSurface", LEAPscintillator_bar, air_x1, OpSciAirSurface);
    new G4LogicalBorderSurface("Air-SciSurface", air_x1, LEAPscintillator_bar, OpAirSciSurface);
    new G4LogicalBorderSurface("Sci-AirSurface", LEAPscintillator_bar, air_x2, OpSciAirSurface);
    new G4LogicalBorderSurface("Air-SciSurface", air_x2, LEAPscintillator_bar, OpAirSciSurface);
    new G4LogicalBorderSurface("Sci-AirSurface", LEAPscintillator_bar, air_y1, OpSciAirSurface);
    new G4LogicalBorderSurface("Air-SciSurface", air_y1, LEAPscintillator_bar, OpAirSciSurface);
    new G4LogicalBorderSurface("Sci-AirSurface", LEAPscintillator_bar, air_y2, OpSciAirSurface);
    new G4LogicalBorderSurface("Air-SciSurface", air_y2, LEAPscintillator_bar, OpAirSciSurface);
    new G4LogicalBorderSurface("Sci-AirSurface", LEAPscintillator_bar, air_top, OpSciAirSurface);
    new G4LogicalBorderSurface("Air-SciSurface", air_top, LEAPscintillator_bar, OpAirSciSurface);
    new G4LogicalBorderSurface("Sci-AirSurface", LEAPscintillator_bar, air_teflon, OpSciAirSurface);
    new G4LogicalBorderSurface("Air-SciSurface", air_teflon, LEAPscintillator_bar, OpAirSciSurface);

    // Optical Pad - Scintillator
    new G4LogicalBorderSurface("Sci-PadSurface", LEAPscintillator_bar, LEAP_Optpad, OpSciPadSurface);
    new G4LogicalBorderSurface("Pad-SciSurface", LEAP_Optpad, LEAPscintillator_bar, OpPadSciSurface);

    // ODM98 adapter
    new G4LogicalBorderSurface("Sci-AirSurface", LEAP_TeflonAdapter, air_teflon, OpTeflonAirSurface);
    new G4LogicalBorderSurface("Air-SciSurface", air_teflon, LEAP_TeflonAdapter, OpAirTeflonSurface);
    if (LEAPScint_conical){
      new G4LogicalBorderSurface("Sci-AirSurface", LEAP_TeflonAdapter, air_teflon_cone, OpTeflonAirSurface);
      new G4LogicalBorderSurface("Air-SciSurface", air_teflon_cone, LEAP_TeflonAdapter, OpAirTeflonSurface);
    }
    new G4LogicalBorderSurface("Sci-PadSurface", LEAP_TeflonAdapter, LEAP_Optpad, OpTeflonAirSurface); // between teflon adapter and opt pad, probably negligible
    new G4LogicalBorderSurface("Pad-SciSurface", LEAP_Optpad, LEAP_TeflonAdapter, OpAirTeflonSurface);
*/
    // add PEEK housing? -> assuming 100% absorbing, does not change the light output result anyway


    // ############################### LEAP with 8mm cylindrical bars ###############################
/*
    G4VPhysicalVolume* LEAPscintillator_bar = new G4PVPlacement(NULL, G4ThreeVector(xy_shift, xy_shift, 0), LEAPCylindricalScintLog_, "LEAPScintillator", DetectorLog_, false, 0);

    G4VPhysicalVolume* LEAP_Optpad = new G4PVPlacement(NULL, G4ThreeVector(xy_shift, xy_shift, -display_split-LEAPScint_length/2.-LEAP_Optpad_thick/2.), LEAPopt_pad_log, "Optpad", DetectorLog_, false, 0);
    G4VPhysicalVolume* LEAP_PMT = new G4PVPlacement(NULL, G4ThreeVector(xy_shift, xy_shift, -3*display_split-LEAPScint_length/2.-LEAP_Optpad_thick-LEAP_PMT_thick/2.), LEAPPMTLog_, "PMT", DetectorLog_, false, 0);
    
    G4VPhysicalVolume* LEAP_TeflonAdapter = new G4PVPlacement(NULL, G4ThreeVector(xy_shift, xy_shift, -2*display_split-LEAPScint_length/2.-LEAPTeflonsquare_thick/2.-0.05), TeflonLog_, "TeflonAdapter", DetectorLog_, false, 0);

    G4VPhysicalVolume* LEAPvikuiti_xy = new G4PVPlacement(NULL, G4ThreeVector(xy_shift, -display_split+xy_shift, 0.05/2.*mm), LEAPCylindricalvikuit_xy_log, "Vikuiti", DetectorLog_, false, 0);
    G4VPhysicalVolume* LEAPvikuiti_top = new G4PVPlacement(NULL, G4ThreeVector(xy_shift, xy_shift, display_split+0.05+(LEAPScint_length+Vikuiti_thickness)/2.*mm), LEAPCylindricalvikuit_z_log, "Vikuiti", DetectorLog_, false, 0);

    G4VPhysicalVolume* air_xy = new G4PVPlacement(NULL, G4ThreeVector(xy_shift, xy_shift, 0.05/2.*mm), LEAPCylindricalair_xy_log, "Air", DetectorLog_, false, 0);
    G4VPhysicalVolume* air_top = new G4PVPlacement(NULL, G4ThreeVector(xy_shift, xy_shift, (LEAPScint_length+0.05)/2.*mm), LEAPCylindricalair_z_log, "Air", DetectorLog_, false, 0);

    // PMT
    new G4LogicalBorderSurface("ResinSiPMSurface", LEAPscintillator_bar, LEAP_PMT, OpResinSiPMSurface);
    new G4LogicalBorderSurface("ResinSiPMSurface", LEAP_PMT, LEAPscintillator_bar, OpSiPMResinSurface);

    // Vikuiti
    new G4LogicalBorderSurface("Vikuiti-AirSurface", LEAPvikuiti_xy, air_xy, OpWrapAirSurface);
    new G4LogicalBorderSurface("Air-VikuitiSurface", air_xy, LEAPvikuiti_xy, OpAirWrapSurface);
    new G4LogicalBorderSurface("Vikuiti-AirSurface", LEAPvikuiti_top, air_top, OpWrapAirSurface);
    new G4LogicalBorderSurface("Air-VikuitiSurface", air_top, LEAPvikuiti_top, OpAirWrapSurface);

    // Scintillator
    new G4LogicalBorderSurface("Sci-AirSurface", LEAPscintillator_bar, air_xy, OpSciAirSurface);
    new G4LogicalBorderSurface("Air-SciSurface", air_xy, LEAPscintillator_bar, OpAirSciSurface);
    new G4LogicalBorderSurface("Sci-AirSurface", LEAPscintillator_bar, air_top, OpSciAirSurface);
    new G4LogicalBorderSurface("Air-SciSurface", air_top, LEAPscintillator_bar, OpAirSciSurface);
    //new G4LogicalBorderSurface("Sci-AirSurface", LEAPscintillator_bar, air_teflon, OpSciAirSurface);
    //new G4LogicalBorderSurface("Air-SciSurface", air_teflon, LEAPscintillator_bar, OpAirSciSurface);

    // Optical Pad - Scintillator
    new G4LogicalBorderSurface("Sci-PadSurface", LEAPscintillator_bar, LEAP_Optpad, OpSciPadSurface);
    new G4LogicalBorderSurface("Pad-SciSurface", LEAP_Optpad, LEAPscintillator_bar, OpPadSciSurface);

    // ODM98 adapter
    //new G4LogicalBorderSurface("Sci-AirSurface", LEAP_TeflonAdapter, air_teflon, OpTeflonAirSurface);
    //new G4LogicalBorderSurface("Air-SciSurface", air_teflon, LEAP_TeflonAdapter, OpAirTeflonSurface);
    new G4LogicalBorderSurface("Sci-PadSurface", LEAP_TeflonAdapter, LEAP_Optpad, OpTeflonAirSurface); // between teflon adapter and opt pad, probably negligible
    new G4LogicalBorderSurface("Pad-SciSurface", LEAP_Optpad, LEAP_TeflonAdapter, OpAirTeflonSurface);
*/



    // ############################### Test for Federco/Neutrino experiment ###############################

// determining the light output for an EJ-200 unwrapped scintillator with dimensions of 3x3x120mm^2

    neutrino_scintillator_bar = new G4PVPlacement(NULL, G4ThreeVector(0, 0, 0), NeutrinoScintLog_, "NeutrinoScintillator", DetectorLog_, false, 0);



  // do this properly by adding the 5 layers of air on 5 faces
    new G4LogicalBorderSurface("Sci-AirSurface", neutrino_scintillator_bar, DetectorLog_, OpSciAirSurface);
    new G4LogicalBorderSurface("Air-SciSurface", DetectorLog_, neutrino_scintillator_bar, OpAirSciSurface);

  // add a perfect detector for readout on one extremity

  // gps: inject light with cos distribution at 20mm from the detector edge, l.o. expected to be 1/3

}
 


