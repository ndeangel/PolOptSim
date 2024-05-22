#ifndef POLAR2Hit_hh
#define POLAR2Hit_hh 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "tls.hh"
#include "G4ThreeVector.hh"

class POLAR2Hit: public G4VHit {
public:
    POLAR2Hit();
    ~POLAR2Hit();

    inline void* operator new(size_t);
    inline void  operator delete(void* aHit);

public:
    // hit information
    G4int             TrackID;
    G4String          ParticleName;
    G4int             ParticleCode;
    G4double          GlobalTime;
    G4int             BarID;
    G4int             ModID;
    G4ThreeVector     LocalPosition;
    G4bool            IsEntering;
    G4bool            IsLeaving;
    G4double          PreMomTheta;
    G4double          PreMomPhi;
    G4double          PostMomTheta;
    G4double          PostMomPhi;
    G4double          PreKinEnergy;
    G4double          PostKinEnergy;
    G4String          ProcessName;
    G4double          EnergyDep;
    G4double          EnergyVis;
    G4double          DeltaTime;

private:
    static G4ThreadLocal G4Allocator<POLAR2Hit>* POLAR2HitAllocator_;

};

typedef G4THitsCollection<POLAR2Hit> POLAR2HitsCollection;

inline void* POLAR2Hit::operator new(size_t) {
    if (POLAR2HitAllocator_ == NULL)
        POLAR2HitAllocator_ = new G4Allocator<POLAR2Hit>();
    return static_cast<void*>(POLAR2HitAllocator_->MallocSingle());
}

inline void POLAR2Hit::operator delete(void* aHit) {
    POLAR2HitAllocator_->FreeSingle(static_cast<POLAR2Hit*>(aHit));
}

#endif
