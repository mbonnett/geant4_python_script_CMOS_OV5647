// EnergySD.cc
#include "EnergySD.hh"

#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"

#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4ios.hh"


EnergySD::EnergySD(const G4String& name
                  //,const G4String& hitsCollectionName
                )
    : G4VSensitiveDetector(name),
      fHitsCollectionID(-1),
      fHitsCollection(NULL)
     {
       collectionName.insert("pixelSDHitsCollection");
     }

EnergySD::~EnergySD() {}

void EnergySD::Initialize(G4HCofThisEvent* HCE) {

    //fHitsCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID(GetName() + "/HitsCollection");
    //HCE->AddHitsCollection(fHitsCollectionID, new G4THitsCollection<EnergyHit>(GetName()));

    // Create hits collection
    //fHitsCollection = new pixelSDHitsCollection(SensitiveDetectorName, collectionName[0]);
    fHitsCollection = new pixelSDHitsCollection(GetName(), collectionName[0]);
    // Add this collection in HCE
    static G4int HCID = -1;
    if (HCID<0) {
      HCID = GetCollectionID(0);
    }
    //G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
    HCE->AddHitsCollection( HCID, fHitsCollection );

}

G4bool EnergySD::ProcessHits(G4Step* step,
                             G4TouchableHistory* /*history*/)
{
    const G4Track* track = step->GetTrack();
    if (track->GetParentID() >= 0) {
    G4double depositedEnergy = step->GetTotalEnergyDeposit();
    G4double emittedEnergy = step->GetPreStepPoint()->GetKineticEnergy();

    EnergyHit* newHit = new EnergyHit();
    newHit->SetDepositedEnergy(depositedEnergy);
    newHit->SetEmittedEnergy(emittedEnergy);

    //G4HCofThisEvent* HCE = step->GetEvent()->GetHCofThisEvent();
    //if (!HCE) return false;

  //  EnergyHit* hitsCollection = static_cast<EnergyHit*>(HCE->GetHC(fHitsCollectionID));
    fHitsCollection->insert(newHit);
    }
    return true;
}
