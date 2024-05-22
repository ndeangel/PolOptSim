#include "POLAR2Run.hh"
#include "POLAR2Hit.hh"
#include "G4SDManager.hh"
#include "POLAR2GlobalConfig.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4RunManager.hh"
#include "g4root.hh"

POLAR2Run::POLAR2Run() {
    totalHits_ = 0;
    totalEvents_ = 0;
}

POLAR2Run::~POLAR2Run() {

}

void POLAR2Run::RecordEvent(const G4Event* anEvent) {
    G4PrimaryVertex* primaryVertex = anEvent->GetPrimaryVertex(0);
    if (primaryVertex == NULL)
        return;

    G4Run::RecordEvent(anEvent);

    POLAR2GlobalConfig* fPOLAR2GlobalConfig = POLAR2GlobalConfig::Instance();

    // calculate and record value
    G4int theEventID = anEvent->GetEventID();
    POLAR2HitsCollection* fPOLAR2HitsCollection = NULL;
    G4int numberOfHits = 0;
    if (!fPOLAR2GlobalConfig->primary_only) {
        G4int collectionID = G4SDManager::GetSDMpointer()->GetCollectionID("POLAR2HitsCollection");
        fPOLAR2HitsCollection = static_cast<POLAR2HitsCollection*>(anEvent->GetHCofThisEvent()->GetHC(collectionID));
        if (fPOLAR2HitsCollection == NULL) {
            G4ExceptionDescription msg;
            msg << "Cannot access hitsCollection: POLAR2HitsCollection with ID " << collectionID;
            G4Exception("POLAR2EventAction::EndOfEventAction", "MyCode0001", FatalException, msg);
        }
        numberOfHits = fPOLAR2HitsCollection->entries();
        if (numberOfHits < 1) return;
        double energy_dep;
        energy_dep = 0.0;

        for (int i = 0; i < numberOfHits; i++) {
            POLAR2Hit* aHit = static_cast<POLAR2Hit*>(fPOLAR2HitsCollection->GetHit(i));
            energy_dep += aHit->EnergyDep;
        }
    }

    G4AnalysisManager* fAnalysisManager = G4AnalysisManager::Instance();
    G4int firstNtupleId = 0;
    if (G4RunManager::GetRunManager()->GetRunManagerType() == G4RunManager::RMType::sequentialRM) {
        firstNtupleId = 1;
    }

    // read primary
    if (anEvent->GetNumberOfPrimaryVertex() > 1) {
        G4ExceptionDescription msg;
        msg << "More than one primary particle for one event.";
        G4Exception("POLAR2EventAction::EndOfEventAction", "MyCode0002", FatalException, msg);
    }
    G4PrimaryParticle* primaryParticle = primaryVertex->GetPrimary(0);
    G4String        pname_primary   = primaryParticle->GetParticleDefinition()->GetParticleName();
    G4ThreeVector   xyz_primary     = primaryVertex->GetPosition();
    G4double        tim_primary     = primaryVertex->GetT0();
    G4double        ene_primary     = primaryParticle->GetKineticEnergy();
    G4ThreeVector   pol_primary     = primaryParticle->GetPolarization();
    G4ThreeVector   direction       = primaryParticle->GetMomentumDirection();
    G4double        theta_primary   = (-direction).theta();
    G4double        phi_primary     = (-direction).phi();
    if (phi_primary < 0)
        phi_primary += twopi;

    // fill ntuple
    // necessary
    fAnalysisManager->FillNtupleIColumn(firstNtupleId + 0, 0 , theEventID               );
    fAnalysisManager->FillNtupleIColumn(firstNtupleId + 0, 1 , numberOfHits             );
    // optional
    if (fPOLAR2GlobalConfig->simout_more) {
        fAnalysisManager->FillNtupleSColumn(firstNtupleId + 0, 2 , pname_primary        );
        fAnalysisManager->FillNtupleFColumn(firstNtupleId + 0, 3 , xyz_primary.x() / cm );
        fAnalysisManager->FillNtupleFColumn(firstNtupleId + 0, 4 , xyz_primary.y() / cm );
        fAnalysisManager->FillNtupleFColumn(firstNtupleId + 0, 5 , xyz_primary.z() / cm );
        fAnalysisManager->FillNtupleDColumn(firstNtupleId + 0, 6 , tim_primary / second );
        fAnalysisManager->FillNtupleFColumn(firstNtupleId + 0, 7 , ene_primary / keV    );
        fAnalysisManager->FillNtupleFColumn(firstNtupleId + 0, 8 , pol_primary.x()      );
        fAnalysisManager->FillNtupleFColumn(firstNtupleId + 0, 9 , pol_primary.y()      );
        fAnalysisManager->FillNtupleFColumn(firstNtupleId + 0, 10, pol_primary.z()      );
        fAnalysisManager->FillNtupleFColumn(firstNtupleId + 0, 11, theta_primary / deg  );
        fAnalysisManager->FillNtupleFColumn(firstNtupleId + 0, 12, phi_primary / deg    );
    }
    fAnalysisManager->AddNtupleRow(firstNtupleId + 0);
    totalEvents_++;

    if (fPOLAR2GlobalConfig->primary_only) return;

    if (fPOLAR2GlobalConfig->event_verbose > 1)
        G4cout << "---> (Record Event) Begin of event: " << theEventID << G4endl;

    // save hits
    for (int i = 0; i < numberOfHits; i++) {
        POLAR2Hit* aHit = static_cast<POLAR2Hit*>(fPOLAR2HitsCollection->GetHit(i));
        // necessary
        fAnalysisManager->FillNtupleIColumn(firstNtupleId + 1, 0 , theEventID                       );
        fAnalysisManager->FillNtupleIColumn(firstNtupleId + 1, 1 , aHit->TrackID                    );
        fAnalysisManager->FillNtupleIColumn(firstNtupleId + 1, 2 , aHit->ParticleCode               );
        fAnalysisManager->FillNtupleIColumn(firstNtupleId + 1, 3 , aHit->BarID                      );
        fAnalysisManager->FillNtupleIColumn(firstNtupleId + 1, 4 , aHit->ModID                      );
        fAnalysisManager->FillNtupleFColumn(firstNtupleId + 1, 5 , aHit->LocalPosition.x() / mm     );
        fAnalysisManager->FillNtupleFColumn(firstNtupleId + 1, 6 , aHit->LocalPosition.y() / mm     );
        fAnalysisManager->FillNtupleFColumn(firstNtupleId + 1, 7 , aHit->LocalPosition.z() / mm     );
        fAnalysisManager->FillNtupleFColumn(firstNtupleId + 1, 8 , aHit->EnergyVis / keV            );
        // optional
        if (fPOLAR2GlobalConfig->simout_more) {
            fAnalysisManager->FillNtupleSColumn(firstNtupleId + 1, 9 , aHit->ParticleName           );
            fAnalysisManager->FillNtupleDColumn(firstNtupleId + 1, 10, aHit->GlobalTime / second    );
            G4int stepStatus = 0;
            if (aHit->IsEntering) stepStatus += 1;
            if (aHit->IsLeaving) stepStatus += 2;
            fAnalysisManager->FillNtupleIColumn(firstNtupleId + 1, 11, stepStatus                   );
            fAnalysisManager->FillNtupleFColumn(firstNtupleId + 1, 12, aHit->PreMomTheta / deg      );
            fAnalysisManager->FillNtupleFColumn(firstNtupleId + 1, 13, aHit->PreMomPhi / deg        );
            fAnalysisManager->FillNtupleFColumn(firstNtupleId + 1, 14, aHit->PostMomTheta / deg     );
            fAnalysisManager->FillNtupleFColumn(firstNtupleId + 1, 15, aHit->PostMomPhi / deg       );
            fAnalysisManager->FillNtupleFColumn(firstNtupleId + 1, 16, aHit->PreKinEnergy / keV     );
            fAnalysisManager->FillNtupleFColumn(firstNtupleId + 1, 17, aHit->PostKinEnergy / keV    );
            fAnalysisManager->FillNtupleSColumn(firstNtupleId + 1, 18, aHit->ProcessName            );
            fAnalysisManager->FillNtupleFColumn(firstNtupleId + 1, 19, aHit->EnergyDep / keV        );
            fAnalysisManager->FillNtupleDColumn(firstNtupleId + 1, 20, aHit->DeltaTime / second     );
        }
        fAnalysisManager->AddNtupleRow(firstNtupleId + 1);
        totalHits_++;
    }

    if (fPOLAR2GlobalConfig->event_verbose > 1)
        G4cout << "---> (Record Event) End of event: " << theEventID << G4endl;
}

void POLAR2Run::Merge(const G4Run* aRun) {
    G4Run::Merge(aRun);

    const POLAR2Run* localRun = static_cast<const POLAR2Run*>(aRun);
    totalHits_ += localRun->totalHits_;
    totalEvents_ += localRun->totalEvents_;

}
