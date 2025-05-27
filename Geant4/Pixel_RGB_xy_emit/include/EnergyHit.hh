// EnergyHit.hh
#include "G4VHit.hh"
#include "G4THitsCollection.hh"

class EnergyHit : public G4VHit {
public:
    EnergyHit();
    virtual ~EnergyHit();

    void SetDepositedEnergy(G4double energy);
    void SetEmittedEnergy(G4double energy);

    G4double GetDepositedEnergy() const;
    G4double GetEmittedEnergy() const;

private:
    G4double fDepositedEnergy;
    G4double fEmittedEnergy;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........B2TrackeroooOO0OOooo......

typedef G4THitsCollection<EnergyHit> pixelSDHitsCollection;
//typedef G4THitsCollection<EnergyHit> GetName();
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//#endif
