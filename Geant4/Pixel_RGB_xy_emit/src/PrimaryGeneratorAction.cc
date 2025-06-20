//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
//#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ChargedGeantino.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include "G4GeneralParticleSource.hh"
//#include "G4Geantino.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
    : G4VUserPrimaryGeneratorAction(),
      fParticleGun(0)
{
    //G4int n_particle = 1;
    //fParticleGun  = new G4ParticleGun(n_particle);
    fParticleGun  = new G4GeneralParticleSource();

    // default particle kinematic


    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* particle = particleTable->FindParticle("chargedgeantino");
    fParticleGun->SetParticleDefinition(particle);
    //fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,0.));
    //fParticleGun->SetParticleEnergy(1*eV);
    //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(1.,0.,0.));

//	G4ParticleDefinition* particleDefinition
//                    = G4ParticleTable::GetParticleTable()->FindParticle("e-");
//	fParticleSource->SetParticleDefinition(particleDefinition);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
/*    G4ParticleDefinition* particle = fParticleGun->GetParticleDefinition();
    if (particle == G4ChargedGeantino::ChargedGeantino())
    {
        //ion
        G4int Z = 38, A = 90;
        G4double ionCharge   = 0.*eplus;
        G4double excitEnergy = 0.*keV;

        G4ParticleDefinition* ion
            = G4IonTable::GetIonTable()->GetIon(Z,A,excitEnergy);
        fParticleGun->SetParticleDefinition(ion);
        fParticleGun->SetParticleCharge(ionCharge);
    }
    */
    /*
      // randomized position
      //
      G4double x0  = 0*um, y0  = 0*um, z0  = 0*um;
      G4double dx0 = 000.*um, dy0 = 000.*um, dz0 = 000.*um;
      //G4double x0  = 4*cm, y0  = 4*cm, z0  = 4*cm;
      //G4double dx0 = 1*cm, dy0 = 1*cm, dz0 = 1*cm;
      x0 += dx0*(G4UniformRand()-0.5);
      y0 += dy0*(G4UniformRand()-0.5);
      z0 += dz0*(G4UniformRand()-0.5);
      fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
     */
    //create vertex
    //
    fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
