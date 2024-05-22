#ifndef POLAR2RunAction_h
#define POLAR2RunAction_h

#include "globals.hh"
#include "G4UserRunAction.hh"

class POLAR2RunAction: public G4UserRunAction {
private:
    G4String output_file_;
    G4bool   fixed_output_;
    static G4int randomSeed_;
    static G4int runId_;

public:
    POLAR2RunAction(G4String the_output_file_, G4bool fixed_name = false);
    ~POLAR2RunAction();

    G4Run* GenerateRun();
    void BeginOfRunAction(const G4Run* aRun);
    void EndOfRunAction(const G4Run* aRun);

};

#endif /* POLAR2RunAction_h */
