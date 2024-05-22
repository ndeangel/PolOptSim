#include "POLAR2EventAction.hh"
#include "POLAR2GlobalConfig.hh"

POLAR2EventAction::POLAR2EventAction() {

}

POLAR2EventAction::~POLAR2EventAction() {

}

void POLAR2EventAction::BeginOfEventAction(const G4Event* anEvent) {
    POLAR2GlobalConfig* fPOLAR2GlobalConfig = POLAR2GlobalConfig::Instance();
    if (fPOLAR2GlobalConfig->event_verbose < 1) return;
    G4int nEvent = anEvent -> GetEventID();
    if((nEvent % 100 == 0) && (nEvent != 0)) {
        G4cout << "INFORMATION: " << nEvent << " event in progress..." << G4endl;
    }
}

void POLAR2EventAction::EndOfEventAction(const G4Event* anEvent) {
    POLAR2GlobalConfig* fPOLAR2GlobalConfig = POLAR2GlobalConfig::Instance();
    if (fPOLAR2GlobalConfig->event_verbose < 1) return;
    G4int nmbEvents = (anEvent -> GetEventID()) + 1;
    if(nmbEvents % 100 == 0)
    {
        G4cout << "INFORMATION: " << nmbEvents << " events processed." << G4endl;
    }
}
