#ifndef POLAR2RUN_HH
#define POLAR2RUN_HH

#include "G4Run.hh"
#include "G4Event.hh"
#include "globals.hh"

class POLAR2Run: public G4Run {
public:
    POLAR2Run();
    ~POLAR2Run();

    void RecordEvent(const G4Event* anEvent);
    void Merge(const G4Run* aRun);

private:
    G4int totalHits_;
    G4int totalEvents_;

public:
    G4int GetTotalHits() const {return totalHits_;}
    G4int GetTotalEvents() const {return totalEvents_;}

};

#endif
