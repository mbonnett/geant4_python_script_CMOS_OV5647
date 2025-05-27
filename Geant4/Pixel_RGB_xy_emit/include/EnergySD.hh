// EnergySD.hh
#ifndef EnergySD_h
#define EnergySD_h 1

#include "G4VSensitiveDetector.hh"
#include "EnergyHit.hh"
///#include "G4HCofThisEvent.hh"
///#include "G4Step.hh"

#include <vector>

class G4Step;
class G4HCofThisEvent;

class EnergySD : public G4VSensitiveDetector {
public:
    EnergySD(const G4String& name
                          //,const G4String& hitsCollectionName
                        );
    virtual ~EnergySD();

    virtual void Initialize(G4HCofThisEvent* HCE);
    virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
  //  virtual void EndOfEvent(G4HCofThisEvent* hitCollection);

private:
    G4int fHitsCollectionID;
    pixelSDHitsCollection* fHitsCollection;
};

#endif
