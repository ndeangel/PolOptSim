#ifndef POLAR2SteppingAction_h
#define POLAR2SteppingAction_h 1

#include "G4UserSteppingAction.hh"

class POLAR2SteppingAction: public G4UserSteppingAction {
public:
    POLAR2SteppingAction();
    ~POLAR2SteppingAction();

    void UserSteppingAction(const G4Step* aStep);

};

#endif
