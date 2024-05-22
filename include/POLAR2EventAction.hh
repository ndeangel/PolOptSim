#ifndef POLAR2EventAction_hh
#define POLAR2EventAction_hh 1

#include "G4UserEventAction.hh"
#include "G4Event.hh"

class POLAR2EventAction: public G4UserEventAction {
public:
    POLAR2EventAction();
    ~POLAR2EventAction();

    void BeginOfEventAction(const G4Event* anEvent);
    void EndOfEventAction(const G4Event* anEvent);

};


#endif
