// EnergyHit.cc
#include "EnergyHit.hh"

EnergyHit::EnergyHit() : G4VHit(), fDepositedEnergy(0.0), fEmittedEnergy(0.0) {}

EnergyHit::~EnergyHit() {}

void EnergyHit::SetDepositedEnergy(G4double energy) {
    fDepositedEnergy = energy;
}

void EnergyHit::SetEmittedEnergy(G4double energy) {
    fEmittedEnergy = energy;
}

G4double EnergyHit::GetDepositedEnergy() const {
    return fDepositedEnergy;
}

G4double EnergyHit::GetEmittedEnergy() const {
    return fEmittedEnergy;
}
