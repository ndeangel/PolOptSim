#ifndef POLAR2SensitiveDetector_hh
#define POLAR2SensitiveDetector_hh 1

#include "G4VSensitiveDetector.hh"
#include "G4EmSaturation.hh"
#include "POLAR2Hit.hh"

class POLAR2SensitiveDetector: public G4VSensitiveDetector {
public:
    POLAR2SensitiveDetector(G4String SDName, G4String HCName);
    ~POLAR2SensitiveDetector();

    void Initialize(G4HCofThisEvent* hitCollection);
    G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);
    void EndOfEvent(G4HCofThisEvent* hitCollection);

private:
    POLAR2HitsCollection * fPOLAR2HitsCollection_;
    G4EmSaturation* fEmSaturation_;

};

#endif
