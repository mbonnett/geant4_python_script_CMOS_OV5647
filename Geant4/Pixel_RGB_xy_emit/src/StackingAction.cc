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
/// \file StackingAction.cc
/// \brief Implementation of the StackingAction class

#include "StackingAction.hh"
#include "RunAction.hh"
//#include "Analysis.hh"

#include "G4Track.hh"
#include "G4NeutrinoE.hh"

#include "G4RunManager.hh"
#include "G4HadronicProcessType.hh"

#include "G4SystemOfUnits.hh"
//#include "G4VProcess.hh"
//#include "G4AnalysisManager.hh"

#include "G4Gamma.hh"

#include <iostream>
#include <fstream>
using namespace std;

#include "EventAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//StackingAction::StackingAction()
//{ }

StackingAction::StackingAction(RunAction* runAction)
    : G4UserStackingAction(),
      fRunAction(runAction)
{ }
//StackingAction::StackingAction(EventAction* eventAction)
  //  : G4UserStackingAction(),
    //  fEventAction(eventAction)
//{ }



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StackingAction::~StackingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ClassificationOfNewTrack
StackingAction::ClassifyNewTrack(const G4Track* track)
{

G4int feventID = (G4EventManager::GetEventManager())->GetConstCurrentEvent()->GetEventID();
    G4String fSource = fRunAction->GetSource();
    G4String fileName = fRunAction->GetFilename();

    fileName += "_";
    fileName += fRunAction->GetSource();



    G4String  emit_allRDecayProducts;
    emit_allRDecayProducts += fileName;
    emit_allRDecayProducts +="_emit_allRDecayProducts";
    emit_allRDecayProducts +=".csv";
    ofstream file_emit_allRDecayProducts (emit_allRDecayProducts, ios::out | ios::app);

    G4String  emit_ParticlesCharge;
    emit_ParticlesCharge += fileName;
    emit_ParticlesCharge +="_emit_ParticlesCharge";
    emit_ParticlesCharge +=".csv";
    ofstream file_emit_ParticlesCharge (emit_ParticlesCharge, ios::out | ios::app);

    G4String  emit_electron;
    emit_electron += fileName;
    emit_electron +="_emit_electron";
    emit_electron +=".csv";
    ofstream file_emit_electron (emit_electron, ios::out | ios::app);

    G4String  emit_positron;
    emit_positron += fileName;
    emit_positron +="_emit_positron";
    emit_positron +=".csv";
    ofstream file_emit_positron (emit_positron, ios::out | ios::app);

    G4String  emit_alpha;
    emit_alpha += fileName;
    emit_alpha +="_emit_alpha";
    emit_alpha +=".csv";
    ofstream file_emit_alpha (emit_alpha, ios::out | ios::app);

    G4String  emit_gamma;
    emit_gamma += fileName;
    emit_gamma +="_emit_gamma";
    emit_gamma +=".csv";
    ofstream file_emit_gamma (emit_gamma, ios::out | ios::app);

    G4String  emit_anti_nu_elec;
    emit_anti_nu_elec += fileName;
    emit_anti_nu_elec +="_emit_anti_nu_elec";
    emit_anti_nu_elec +=".csv";
    ofstream file_emit_anti_nu_elec (emit_anti_nu_elec, ios::out | ios::app);


//if (track->GetTrackID() != 1)





    //keep primary particle
    if (track->GetParentID() == 0) return fUrgent;

//////////////////////////////////////////////////////////////
    //G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

// Run* run = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
// run->CountParticles(track->GetDefinition());

    //secondary particles only
    //if (track->GetParentID() > 0) return;
    //if (track->GetTrackID() == 1) return;

    const G4ParticleDefinition* particle = track->GetParticleDefinition();
    G4String name   = particle->GetParticleName();
    //G4int pid       = particle->GetPDGEncoding();
    //G4int Z         = particle->GetAtomicNumber();
    //G4int A         = particle->GetAtomicMass();
    G4double charge = particle->GetPDGCharge();
    G4double energy = track->GetKineticEnergy();
    //G4double time   = track->GetGlobalTime();
  //  G4double weight = track->GetWeight();
    G4double totalenergy = track->GetTotalEnergy () ;
    G4ThreeVector pos_xyz = track->GetPosition ();
    //  G4int N_step  = track->       //GetTrackID(); ///GetParentID();
    //G4double PosX = pos_xyz.x();
    //G4double PosY = pos_xyz.y();
    //G4double PosZ = pos_xyz.z();

    //G4cout << "ParticleName: " <<name << " PDGEncoding: "<<pid<<G4endl;

    /*  G4cout << "Track # "
             << track->GetTrackID() << " of " << name
             << " E(MeV)= " << energy/MeV
             << " produced by Track ID= " << track->GetParentID()
             << G4endl;

      if (track ->GetTrackID() != 1 && track -> GetParentID() == 1) {
          //std::cout << "Particle is a secondary " <<name << " alltotalenergy: " <<totalenergy << " energy: "<<energy<< std::endl;
          //if ((name == "e-"))
         G4cout << "Particle is a secondary but parent was a primary " <<name << " energy: "<<energy<< G4endl;
      }
    */



//if (track->GetParentID() > 0){
//G4cout << "Particle is a secondary: " << name<<G4endl;
//}
//photon (optical) G4OpticalPhoton opticalphoton (0)

    G4int processType = track->GetCreatorProcess()->GetProcessSubType();

    if(fSource=="e" || fSource=="alpha") {

//G4cout << "ParticleName: " <<name << " PDGEncoding: "<<pid<<G4endl;
//if (charge < 3.) {
G4int id ;
if (charge != 0.)
{
    //fill ntuple id = 0
    id = 3;
    //analysisManager->FillNtupleDColumn(id,0, double(pid));
    //analysisManager->FillNtupleDColumn(id,1, energy);
    //analysisManager->FillNtupleDColumn(id,2, time/s);
    //analysisManager->AddNtupleRow(id);
    if (file_emit_ParticlesCharge.is_open()) {
        file_emit_ParticlesCharge << energy << "\t" <<totalenergy <<"\t"<<
feventID<<endl;
    }
    file_emit_ParticlesCharge.close();

    //analysisManager->FillH1(2, energy);
    // //analysisManager->SetH1Ascii(2, true);
}

if (name == "e-")
{
    id = 4;
    //analysisManager->FillNtupleDColumn(id,0, energy);
    //analysisManager->FillNtupleDColumn(id,1, time/s);
    //analysisManager->AddNtupleRow(id);
    // if (file_spectrm.is_open()){file_spectrm << energy << "\t" <<totalenergy << endl;}
    //  G4cout << "e-totalenergy: " << "\t" <<totalenergy << " energy: "<<energy<<G4endl;

    if (file_emit_electron.is_open()) {
        file_emit_electron << energy << "\t" <<totalenergy << "\t"<<
feventID<<endl;
    }
    file_emit_electron.close();

    //analysisManager->FillH1(3, energy);
    // //analysisManager->SetH1Ascii(3, true);
}

if (name == "e+")
{
    id = 5;
    //analysisManager->FillNtupleDColumn(id,0, energy);
    //analysisManager->FillNtupleDColumn(id,1, time/s);
    //analysisManager->AddNtupleRow(id);
    if (file_emit_positron.is_open()) {
        file_emit_positron << energy << "\t" <<totalenergy << endl;
    }
    file_emit_positron.close();

    ////analysisManager->FillH1(3, energy, weight);
}

if (name == "alpha")
{
    id = 6;
    //analysisManager->FillNtupleDColumn(id,0, energy);
    //analysisManager->FillNtupleDColumn(id,1, time/s);
    //analysisManager->AddNtupleRow(id);
    if (file_emit_alpha.is_open()) {
        file_emit_alpha << energy << "\t" <<totalenergy << "\t"<<
feventID<<endl;
    }
    file_emit_alpha.close();

    ////analysisManager->FillH1(3, energy, weight);
}

if (name == "gamma")
{
    id = 7;
    //analysisManager->FillNtupleDColumn(id,0, energy);
    //analysisManager->FillNtupleDColumn(id,1, time/s);
    //analysisManager->AddNtupleRow(id);
    if (file_emit_gamma.is_open()) {
        file_emit_gamma << energy << "\t" <<totalenergy << "\t"<<
feventID<<endl;
    }
    file_emit_gamma.close();

    ////analysisManager->FillH1(3, energy, weight);
}

    }

    else{
    if (processType == fRadioactiveDecay)
    {
        G4int id = 2;
        //analysisManager->FillNtupleDColumn(id,0, double(pid));
        //analysisManager->FillNtupleDColumn(id,1, double(Z));
        //analysisManager->FillNtupleDColumn(id,2, double(A));
        //analysisManager->FillNtupleDColumn(id,3, energy);
        //analysisManager->FillNtupleDColumn(id,4, time/s);
        //analysisManager->FillNtupleDColumn(id,5, weight);
        //analysisManager->FillNtupleDColumn(id,6, PosX);
        //analysisManager->FillNtupleDColumn(id,7, PosY);
        //analysisManager->FillNtupleDColumn(id,8, PosZ);
        //analysisManager->AddNtupleRow(id);
        //  G4cout <<name << " alltotalenergy: " <<totalenergy << " energy: "<<energy<<G4endl;
        if (file_emit_allRDecayProducts.is_open()) {
            file_emit_allRDecayProducts << energy << "\t" <<totalenergy << "\t"<<
feventID<<endl;
        }
        file_emit_allRDecayProducts.close();

        //if (charge < 3.) {
        if (charge != 0.)
        {
            //fill ntuple id = 0
            id = 3;
            //analysisManager->FillNtupleDColumn(id,0, double(pid));
            //analysisManager->FillNtupleDColumn(id,1, energy);
            //analysisManager->FillNtupleDColumn(id,2, time/s);
            //analysisManager->AddNtupleRow(id);
            if (file_emit_ParticlesCharge.is_open()) {
                file_emit_ParticlesCharge << energy << "\t" <<totalenergy << "\t"<<
feventID<<endl;
            }
            file_emit_ParticlesCharge.close();

            //analysisManager->FillH1(2, energy);
            // //analysisManager->SetH1Ascii(2, true);
        }

        if (name == "e-")
        {
            id = 4;
            //analysisManager->FillNtupleDColumn(id,0, energy);
            //analysisManager->FillNtupleDColumn(id,1, time/s);
            //analysisManager->AddNtupleRow(id);
            // if (file_spectrm.is_open()){file_spectrm << energy << "\t" <<totalenergy << endl;}
            //  G4cout << "e-totalenergy: " << "\t" <<totalenergy << " energy: "<<energy<<G4endl;

            if (file_emit_electron.is_open()) {
                file_emit_electron << energy << "\t" <<totalenergy << "\t"<<
feventID<<endl;
            }
            file_emit_electron.close();

            //analysisManager->FillH1(3, energy);
            // //analysisManager->SetH1Ascii(3, true);
        }

        if (name == "e+")
        {
            id = 5;
            //analysisManager->FillNtupleDColumn(id,0, energy);
            //analysisManager->FillNtupleDColumn(id,1, time/s);
            //analysisManager->AddNtupleRow(id);
            if (file_emit_positron.is_open()) {
                file_emit_positron << energy << "\t" <<totalenergy << "\t"<<
feventID<<endl;
            }
            file_emit_positron.close();

            ////analysisManager->FillH1(3, energy, weight);
        }

        if (name == "alpha")
        {
            id = 6;
            //analysisManager->FillNtupleDColumn(id,0, energy);
            //analysisManager->FillNtupleDColumn(id,1, time/s);
            //analysisManager->AddNtupleRow(id);
            if (file_emit_alpha.is_open()) {
                file_emit_alpha << energy << "\t" <<totalenergy << "\t"<<
feventID<<endl;
            }
            file_emit_alpha.close();

            ////analysisManager->FillH1(3, energy, weight);
        }

        if (name == "gamma")
        {
            id = 7;
            //analysisManager->FillNtupleDColumn(id,0, energy);
            //analysisManager->FillNtupleDColumn(id,1, time/s);
            //analysisManager->AddNtupleRow(id);
            if (file_emit_gamma.is_open()) {
                file_emit_gamma << energy << "\t" <<totalenergy << "\t"<<
feventID<<endl;
            }
            file_emit_gamma.close();

            ////analysisManager->FillH1(3, energy, weight);
        }

        if (name == "anti_nu_e")
        {
            id = 8;
            //analysisManager->FillNtupleDColumn(id,0, energy);
            //analysisManager->FillNtupleDColumn(id,1, time/s);
            //analysisManager->AddNtupleRow(id);
            if (file_emit_anti_nu_elec.is_open()) {
                file_emit_anti_nu_elec << energy << "\t" <<totalenergy << "\t"<<
feventID<<endl;
            }
            file_emit_anti_nu_elec.close();

            //analysisManager->FillH1(4, energy);
            ////analysisManager->SetH1Ascii(4, true);
        }

        if(fSource=="Sr90") {
            if (name == "Y90") //sr90
            {
                id = 9;
                //analysisManager->FillNtupleDColumn(id,0, energy);
                //analysisManager->FillNtupleDColumn(id,1, time/s);
                //analysisManager->AddNtupleRow(id);
                G4String  emit_Y90;
                emit_Y90 += fileName;
                emit_Y90 +="_emit_Y90";
                emit_Y90 +=".csv";
                ofstream file_emit_Y90 (emit_Y90, ios::out | ios::app);
                if (file_emit_Y90.is_open()) {
                    file_emit_Y90 << energy << "\t" <<totalenergy << "\t"<<
feventID<<endl;
                }
                file_emit_Y90.close();

                ////analysisManager->FillH1(3, energy, weight);
            }

            if (name == "Zr90") //sr90
            {
                id = 10;
                //analysisManager->FillNtupleDColumn(id,0, energy);
                //analysisManager->FillNtupleDColumn(id,1, time/s);
                //analysisManager->AddNtupleRow(id);
                G4String  emit_Zr90;
                emit_Zr90 += fileName;
                emit_Zr90 +="_emit_Zr90";
                emit_Zr90 +=".csv";
                ofstream file_emit_Zr90 (emit_Zr90, ios::out | ios::app);
                if (file_emit_Zr90.is_open()) {
                    file_emit_Zr90 << energy << "\t" <<totalenergy << "\t"<<
feventID<<endl;
                }
                file_emit_Zr90.close();
                ////analysisManager->FillH1(3, energy, weight);
            }
        }

        if(fSource=="Cs137") {
            if (name == "Ba137[661.659]") //cs137(1000561371)
            {
                id = 9;
                //analysisManager->FillNtupleDColumn(id,0, energy);
                //analysisManager->FillNtupleDColumn(id,1, time/s);
                //analysisManager->AddNtupleRow(id);
                G4String  emit_Ba137_661_659;
                emit_Ba137_661_659 += fileName;
                emit_Ba137_661_659 +="_emit_emit_Ba137_661_659";
                emit_Ba137_661_659 +=".csv";
                ofstream file_emit_Ba137_661_659 (emit_Ba137_661_659, ios::out | ios::app);
                if (file_emit_Ba137_661_659.is_open()) {
                    file_emit_Ba137_661_659 << energy << "\t" <<totalenergy << "\t"<<
feventID<<endl;
                }
                file_emit_Ba137_661_659.close();
                ////analysisManager->FillH1(3, energy, weight);
            }

            if (name == "Ba137") //cs137(1000250550)
            {
                id = 10;
                //analysisManager->FillNtupleDColumn(id,0, energy);
                //analysisManager->FillNtupleDColumn(id,1, time/s);
                //analysisManager->AddNtupleRow(id);
                G4String  emit_Ba137;
                emit_Ba137 += fileName;
                emit_Ba137 +="_emit_Ba137";
                emit_Ba137 +=".csv";
                ofstream file_emit_Ba137 (emit_Ba137, ios::out | ios::app);
                if (file_emit_Ba137.is_open()) {
                    file_emit_Ba137 << energy << "\t" <<totalenergy << "\t"<<
feventID<<endl;
                }
                file_emit_Ba137.close();

                ////analysisManager->FillH1(3, energy, weight);
            }
        }
        if(fSource=="Fe55") {
            if (name == "Mn55") //Fe55(1000561371)
            {
                id = 9;
                //analysisManager->FillNtupleDColumn(id,0, energy);
                //analysisManager->FillNtupleDColumn(id,1, time/s);
                //analysisManager->AddNtupleRow(id);
                G4String  emit_Mn55;
                emit_Mn55 += fileName;
                emit_Mn55 +="_emit_Mn55";
                emit_Mn55 +=".csv";
                ofstream file_emit_Mn55 (emit_Mn55, ios::out | ios::app);
                if (file_emit_Mn55.is_open()) {
                    file_emit_Mn55 << energy << "\t" <<totalenergy << "\t"<<
feventID<<endl;
                }
                file_emit_Mn55.close();
                //analysisManager->FillH1(3, energy, weight);
            }
        }
    }
//    file_spectrm.close();

   }
//////////////////////////////////////////////////////////////////////////////
//if(track -> GetDefinition() == G4Gamma::Definition()) {
    //  G4cout << "Particle is a secondary Photons: " << name<<G4endl;
//}
//////////////////////////////////////////////////////////////////////////////

    //kill secondary neutrino
    if (track->GetDefinition() == G4NeutrinoE::NeutrinoE()) return fKill;
    else return fUrgent;
//return fUrgent;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
