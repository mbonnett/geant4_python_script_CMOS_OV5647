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
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"

#include "G4Step.hh"
//#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "RunAction.hh"

//#include <G4Electron.hh>
#include <G4SystemOfUnits.hh>
///#include "Analysis.hh"
//#include "G4AnalysisManager.hh"

#include <iostream>
#include <fstream>
using namespace std;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(EventAction* eventAction)
//SteppingAction::SteppingAction(RunAction* runAction)
    : G4UserSteppingAction(),
      //fRunAction(runAction)
      fEventAction(eventAction)
      //fScoringVolume(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{

  G4int feventID = fEventAction->feventID;
  //G4int fPlx_x = fEventAction->fPlx_x;
  //G4int fPlx_y = fEventAction->fPlx_y;
  //G4double fEdep_plx = fEventAction->fEdep_plx;

    // Collect energy and track length step by step
    G4String procName = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
    G4String particleName = step->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName();
    G4int parent_id = step->GetTrack()->GetParentID();
    G4int track_id = step->GetTrack()->GetTrackID();

    //G4StepPoint* prepoint = step->GetPreStepPoint();
    //G4StepPoint* postpoint = step->GetPostStepPoint();

    //G4ThreeVector pos1 = prepoint->GetPosition();
    //  G4ThreeVector pos2 = postpoint->GetPosition();
    G4double posX = step->GetPreStepPoint()->GetPosition().x()/um;
    G4double posY = step->GetPreStepPoint()->GetPosition().y()/um;
    G4double posZ = step->GetPreStepPoint()->GetPosition().z()/um;
    //G4int step_id = step->GetTrack()->GetCurrentStepNumber();


    // collect energy deposited in this step
    G4double edep = step->GetTotalEnergyDeposit();
    G4double ekinpre = step->GetPreStepPoint()->GetKineticEnergy();
    G4double ekinpost = step->GetPostStepPoint()->GetKineticEnergy();
    G4double deltaEnergy = step->GetDeltaEnergy ()/keV;
    G4double ekin = step->GetTrack()->GetKineticEnergy();

    //G4double trackLength = step->GetTrack()->GetTrackLength()/um;
  // step length
    G4double stepLength = step->GetStepLength()/um;
      //G4double stepLength = 0.;
    //if ( step->GetTrack()->GetDefinition()->GetPDGCharge() != 0. ) {
    //  G4double stepLength = step->GetStepLength();   }
    G4double dEdL=abs(deltaEnergy/stepLength);
    //G4double dEdL=((edep/keV)/stepLength);

    // get volume of the current step
    G4VPhysicalVolume* volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
//    const G4ParticleDefinition* particle = step->GetTrack()->GetParticleDefinition();

///////////////////////////////////////////////////////////////////////////////////////////////////////////
    G4String fSource = fEventAction->fSource;


    G4int N_pxl_x = fEventAction->fN_pxl_x;
    G4int N_pxl_y = fEventAction->fN_pxl_y;
    G4double Size_pxl_x = fEventAction->fSize_pxl_x/um;
    G4double Size_pxl_y = fEventAction->fSize_pxl_y/um;
    G4bool SaveProcPart = fEventAction->fSaveProcPart;
//G4cout<<" Process_Name: "<<SaveProcPart<<"    "<< fEventAction->fSaveProcPart<<G4endl;

    G4double length_x = N_pxl_x*Size_pxl_x;
    G4double length_y = N_pxl_y*Size_pxl_y;
    G4double pix_x=(posX+0.5*length_x)/(Size_pxl_x + 0.0000001) ;
    G4double pix_y=(posY+0.5*length_y)/(Size_pxl_y + 0.0000001);
//G4double pix_x=trunc(posX/Size_pxl_x) + (N_pxl_x - 1)*0.5;
//G4double pix_y=trunc(posY/Size_pxl_y) + (N_pxl_y - 1)*0.5;


    G4String fileName = fEventAction->fFilename;
    fileName += "_";
    fileName += fEventAction->fSource;

////////////////////////////////////////////////////////////////////////////////
if(SaveProcPart == 1){
    //G4cout<<" Process_Name: "<<procName<<" Particles_Name: "<<particleName<<G4endl;
    G4String Process_Particles_Name;
    Process_Particles_Name += fileName;
    Process_Particles_Name +="_Process_Particles_Name";
    Process_Particles_Name +=".csv";
    ofstream file_Process_Particles_Name (Process_Particles_Name, ios::out | ios::app);
    if (file_Process_Particles_Name.is_open()) {
        file_Process_Particles_Name<<
                edep/keV<<"\t"<<ekinpost/MeV<<"\t"<<procName<<"\t"<<particleName<<"\t"<<parent_id<<"\t"<<track_id<<"\t"<<
feventID<<endl;
    }
    file_Process_Particles_Name.close();
}
///////??????????????Fuente?????????///////
/*
if (volume && (volume->GetName() != "World") &&
              (volume->GetName() != "Detec_SiO2_layers") &&
              (volume->GetName() != "DetecBase_Si") &&
              (volume->GetName() != "DetecBase_Ss") &&
//(volume->GetName() != "DetecBase") &&
              (volume->GetName() != "Detector") &&

              (volume->GetName() != "Sshield_out") &&
              (volume->GetName() != "cap_Sshield_up") &&

        //  (volume->GetName() != "cap_Sshield_top_up") &&
          //(volume->GetName() != "cap_Stainless_steel") &&

        //  (volume->GetName() != "be_win") &&
          //(volume->GetName() != "cu_shield") &&

              (volume->GetName() != "cap_Sshield_down") ){

    if(particleName == "e-" ) {
        fEventAction->AddE_all_electron_produc(edep, ekinpost);
        //G4cout<<" Process_Name: "<<procName<<"  : "<<particleName<<G4endl;
        G4String All_electron_produc;
        All_electron_produc += fileName;
        All_electron_produc +="_All_electron_produc";
        All_electron_produc +=".csv";
        ofstream file_All_electron_produc (All_electron_produc, ios::out | ios::app);
        if (file_All_electron_produc.is_open()) {
            file_All_electron_produc<<
                    edep/keV<<"\t"<<ekinpost/MeV<<"\t"<<procName<<"\t"<<parent_id<<"\t"<<track_id<<"\t"<<
feventID<<endl;
        }
        file_All_electron_produc.close();
    }

    if(particleName == "gamma" ) {
        //if(particle->GetParticleName() == "gamma" ) {
        fEventAction->AddE_all_gamma_produc(edep, ekinpost);
        //G4cout<<" Process_Name: "<<procName<<" Particles_Name: "<<particleName<<G4endl;
        G4String All_gamma_produc;
        All_gamma_produc += fileName;
        All_gamma_produc +="_All_gamma_produc";
        All_gamma_produc +=".csv";
        ofstream file_All_gamma_produc (All_gamma_produc, ios::out | ios::app);
        if (file_All_gamma_produc.is_open()) {
            file_All_gamma_produc<<
                    edep/keV<<"\t"<<ekinpost/MeV<<"\t"<<procName<<"\t"<<parent_id<<"\t"<<track_id<<"\t"<<
feventID<<endl;
        }
        file_All_gamma_produc.close();
    }
}
*/
///////??????????????Fuente?????????///////

//////////////////////////////////////////////////////////////////////////////////
/////Detector
    if (volume &&  (volume->GetName() == "Detector") ) {

        if(edep>=0. ) {
            //fEventAction->AddE_all_particle_arrived_detector(edep, ekinpost);
            G4String  edep_all_particle_arrived_detector;
            edep_all_particle_arrived_detector += fileName;
            edep_all_particle_arrived_detector +="_edep_all_particle_arrived_detector";
            edep_all_particle_arrived_detector +=".csv";
            ofstream file_edep_all_particle_arrived_detector (edep_all_particle_arrived_detector, ios::out | ios::app);
            if (file_edep_all_particle_arrived_detector.is_open()) {
                file_edep_all_particle_arrived_detector<<
                          trunc((pix_x))<<"\t"<<
                          trunc((pix_y))<<"\t"<<
                          edep/keV<<"\t"<<
                          ekinpost/MeV<<"\t"<<
                          ekinpre<<"\t"<<
                          ekin<<"\t"<<
                          dEdL<<"\t"<<
                          feventID<<"\t"<<
                          particleName<<"\t"<<
                          procName<<"\t"<<
                          parent_id<<"\t"<<
                          track_id<<"\t"<<
                          posX<<"\t"<<
                          posY<<"\t"<<
                          posZ<<"\t"<<
                          endl;

            }
            file_edep_all_particle_arrived_detector.close();
        }

        if(particleName == "e-"  && edep>=0.) {
            fEventAction->AddE_all_electron_arrived_detector(edep, ekinpost, stepLength);///no dexdx
            G4String  edep_all_electron_arrived_detector;
            edep_all_electron_arrived_detector += fileName;
            edep_all_electron_arrived_detector +="_edep_all_electron_arrived_detector";
            edep_all_electron_arrived_detector +=".csv";
            ofstream file_edep_all_electron_arrived_detector (edep_all_electron_arrived_detector, ios::out | ios::app);
            if (file_edep_all_electron_arrived_detector.is_open()) {
                file_edep_all_electron_arrived_detector<<
                          trunc((pix_x))<<"\t"<<
                          trunc((pix_y))<<"\t"<<
                          edep/keV<<"\t"<<
                          ekinpost/MeV<<"\t"<<
                          ekinpre<<"\t"<<
                          ekin<<"\t"<<
                          dEdL<<"\t"<<
                          feventID<<"\t"<<
                          //posX<<"\t"<<posY<<"\t"<<posZ<<"\t"<<
                          //round((pix_x))<<"\t"<<
                          //round((pix_y))<<"\t"<<
                          procName<<"\t"<<
                          parent_id<<"\t"<<
                          track_id<<"\t"<<

                          (posX+0.5*length_x)<<"\t"<<
                          (posY+0.5*length_y)<<"\t"<<
                          posZ<<"\t"<<
                          //fPlx_x<<"\t"<<fPlx_y<<"\t"<<fEdep_plx<<"\t"<<

                          endl;
                //  file_edep_all_electron_arrived_detector<<edep/keV<<"\t"<<ekin/MeV<<"\t"<<procName<<"\t"<<particleName<<"\t"<<ekinpost<<"\t"<<ekinpost<<"\t"<<feventID<<endl;
//file_edep_all_electron_arrived_detector<<edep/keV<<"\t"<<ekin/MeV<<"\t"<<fEventAction->fEdep_plx<<"\t"<<fEventAction->fPlx_x<<"\t"<<fEventAction->fPlx_y<<"\t"<<track_id<<"\t"<<feventID<<endl;
            }
            file_edep_all_electron_arrived_detector.close();
        }
        ///////////////////////////////////////////////
        if(particleName == "anti_nu_e"  && edep>=0.) {
    fEventAction->AddE_all_anti_nu_e_arrived_detector(edep, ekinpost, stepLength);///no dexdx
    G4String  edep_all_anti_nu_e_arrived_detector;
    edep_all_anti_nu_e_arrived_detector += fileName;
    edep_all_anti_nu_e_arrived_detector +="_edep_all_anti_nu_e_arrived_detector";
    edep_all_anti_nu_e_arrived_detector +=".csv";
    ofstream file_edep_all_anti_nu_e_arrived_detector (edep_all_anti_nu_e_arrived_detector, ios::out | ios::app);
    if (file_edep_all_anti_nu_e_arrived_detector.is_open()) {
        file_edep_all_anti_nu_e_arrived_detector<<
                  trunc((pix_x))<<"\t"<<
                  trunc((pix_y))<<"\t"<<
                  edep/keV<<"\t"<<
                  ekinpost/MeV<<"\t"<<
                  ekinpre/MeV<<"\t"<<
                  edep/ekinpre*100<<"\t"<<
                  dEdL<<"\t"<<
                  feventID<<"\t"<<
                  //posX<<"\t"<<posY<<"\t"<<posZ<<"\t"<<
                  //round((pix_x))<<"\t"<<
                  //round((pix_y))<<"\t"<<
                  procName<<"\t"<<
                  parent_id<<"\t"<<
                  track_id<<"\t"<<

                  (posX+0.5*length_x)<<"\t"<<
                  (posY+0.5*length_y)<<"\t"<<
                  posZ<<"\t"<<
                  //fPlx_x<<"\t"<<fPlx_y<<"\t"<<fEdep_plx<<"\t"<<

                  endl;
        //  file_edep_all_anti_nu_e_arrived_detector<<edep/keV<<"\t"<<ekin/MeV<<"\t"<<procName<<"\t"<<particleName<<"\t"<<ekinpost<<"\t"<<ekinpost<<"\t"<<feventID<<endl;
//file_edep_all_anti_nu_e_arrived_detector<<edep/keV<<"\t"<<ekin/MeV<<"\t"<<fEventAction->fEdep_plx<<"\t"<<fEventAction->fPlx_x<<"\t"<<fEventAction->fPlx_y<<"\t"<<track_id<<"\t"<<feventID<<endl;
    }
    file_edep_all_anti_nu_e_arrived_detector.close();
}
        ///////////////////////////////////////////////
        if(particleName == "alpha"  && edep>=0.) {
            fEventAction->AddE_all_alpha_arrived_detector(edep, ekinpost, stepLength);///no dexdx
            G4String  edep_all_alpha_arrived_detector;
            edep_all_alpha_arrived_detector += fileName;
            edep_all_alpha_arrived_detector +="_edep_all_alpha_arrived_detector";
            edep_all_alpha_arrived_detector +=".csv";
            ofstream file_edep_all_alpha_arrived_detector (edep_all_alpha_arrived_detector, ios::out | ios::app);
            if (file_edep_all_alpha_arrived_detector.is_open()) {
                file_edep_all_alpha_arrived_detector<<
                          trunc((pix_x))<<"\t"<<
                          trunc((pix_y))<<"\t"<<
                          edep/keV<<"\t"<<
                          ekinpost/MeV<<"\t"<<
                          ekinpre/MeV<<"\t"<<
                          edep/ekinpre*100<<"\t"<<
                          dEdL<<"\t"<<
                          feventID<<"\t"<<
                          //posX<<"\t"<<posY<<"\t"<<posZ<<"\t"<<
                          //round((pix_x))<<"\t"<<
                          //round((pix_y))<<"\t"<<
                          procName<<"\t"<<
                          parent_id<<"\t"<<
                          track_id<<"\t"<<

                          (posX+0.5*length_x)<<"\t"<<
                          (posY+0.5*length_y)<<"\t"<<
                          posZ<<"\t"<<
                          //fPlx_x<<"\t"<<fPlx_y<<"\t"<<fEdep_plx<<"\t"<<

                          endl;
                //  file_edep_all_alpha_arrived_detector<<edep/keV<<"\t"<<ekin/MeV<<"\t"<<procName<<"\t"<<particleName<<"\t"<<ekinpost<<"\t"<<ekinpost<<"\t"<<feventID<<endl;
//file_edep_all_alpha_arrived_detector<<edep/keV<<"\t"<<ekin/MeV<<"\t"<<fEventAction->fEdep_plx<<"\t"<<fEventAction->fPlx_x<<"\t"<<fEventAction->fPlx_y<<"\t"<<track_id<<"\t"<<feventID<<endl;
            }
            file_edep_all_alpha_arrived_detector.close();
        }
        ///////////////////////////////////////////////
        if(particleName == "gamma"  && edep>=0.) {
            fEventAction->AddE_all_gamma_arrived_detector(edep, ekinpost, stepLength);///no dexdx
            G4String  edep_all_gamma_arrived_detector;
            edep_all_gamma_arrived_detector += fileName;
            edep_all_gamma_arrived_detector +="_edep_all_gamma_arrived_detector";
            edep_all_gamma_arrived_detector +=".csv";
            ofstream file_edep_all_gamma_arrived_detector (edep_all_gamma_arrived_detector, ios::out | ios::app);
            if (file_edep_all_gamma_arrived_detector.is_open()) {
                        file_edep_all_gamma_arrived_detector<<
                        trunc((pix_x))<<"\t"<<
                        trunc((pix_y))<<"\t"<<
                        edep/keV<<"\t"<<
                        ekinpost/MeV<<"\t"<<
                        ekinpre/MeV<<"\t"<<
                        edep/ekinpre*100<<"\t"<<
                        dEdL<<"\t"<<
                        feventID<<"\t"<<
                        //posX<<"\t"<<posY<<"\t"<<posZ<<"\t"<<
                        //round((pix_x))<<"\t"<<
                        //round((pix_y))<<"\t"<<
                        procName<<"\t"<<
                        parent_id<<"\t"<<
                        track_id<<"\t"<<

                        (posX+0.5*length_x)<<"\t"<<
                        (posY+0.5*length_y)<<"\t"<<
                        posZ<<"\t"<<
                        //fPlx_x<<"\t"<<fPlx_y<<"\t"<<fEdep_plx<<"\t"<<

                        endl;
                //  file_edep_all_gamma_arrived_detector<<edep/keV<<"\t"<<ekin/MeV<<"\t"<<procName<<"\t"<<particleName<<"\t"<<ekinpost<<"\t"<<ekinpost<<"\t"<<feventID<<endl;
//file_edep_all_gamma_arrived_detector<<edep/keV<<"\t"<<ekin/MeV<<"\t"<<fEventAction->fEdep_plx<<"\t"<<fEventAction->fPlx_x<<"\t"<<fEventAction->fPlx_y<<"\t"<<track_id<<"\t"<<feventID<<endl;
            }
            file_edep_all_gamma_arrived_detector.close();
        }
        ////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////
        if(edep==0. ) {
            //fEventAction->AddE_zero_all_particle_arrived_detector(edep, ekinpost);
            G4String  edep_zero_all_particle_arrived_detector;
            edep_zero_all_particle_arrived_detector += fileName;
            edep_zero_all_particle_arrived_detector +="_edep_zero_all_particle_arrived_detector";
            edep_zero_all_particle_arrived_detector +=".csv";
            ofstream file_edep_zero_all_particle_arrived_detector (edep_zero_all_particle_arrived_detector, ios::out | ios::app);
            if (file_edep_zero_all_particle_arrived_detector.is_open()) {
                file_edep_zero_all_particle_arrived_detector<<
                          trunc((pix_x))<<"\t"<<
                          trunc((pix_y))<<"\t"<<
                          edep/keV<<"\t"<<
                          ekinpost/MeV<<"\t"<<
                          ekinpre<<"\t"<<
                          ekin<<"\t"<<
                          dEdL<<"\t"<<
                          feventID<<"\t"<<
                          particleName<<"\t"<<
                          procName<<"\t"<<
                          parent_id<<"\t"<<
                          track_id<<"\t"<<
                          posX<<"\t"<<
                          posY<<"\t"<<
                          posZ<<"\t"<<
                          endl;

            }
            file_edep_zero_all_particle_arrived_detector.close();
        }
        //////////////////////////////////////////////
        if(particleName == "e-"  && edep==0.) {
            fEventAction->AddE_zero_all_electron_arrived_detector(edep, ekinpost, stepLength);///no dexdx
            G4String  edep_zero_all_electron_arrived_detector;
            edep_zero_all_electron_arrived_detector += fileName;
            edep_zero_all_electron_arrived_detector +="_edep_zero_all_electron_arrived_detector";
            edep_zero_all_electron_arrived_detector +=".csv";
            ofstream file_edep_zero_all_electron_arrived_detector (edep_zero_all_electron_arrived_detector, ios::out | ios::app);
            if (file_edep_zero_all_electron_arrived_detector.is_open()) {
                file_edep_zero_all_electron_arrived_detector<<
                          trunc((pix_x))<<"\t"<<
                          trunc((pix_y))<<"\t"<<
                          edep/keV<<"\t"<<
                          ekinpost/MeV<<"\t"<<
                          ekinpre<<"\t"<<
                          ekin<<"\t"<<
                          dEdL<<"\t"<<
                          feventID<<"\t"<<
                          //posX<<"\t"<<posY<<"\t"<<posZ<<"\t"<<
                          //round((pix_x))<<"\t"<<
                          //round((pix_y))<<"\t"<<
                          procName<<"\t"<<
                          parent_id<<"\t"<<
                          track_id<<"\t"<<

                          (posX+0.5*length_x)<<"\t"<<
                          (posY+0.5*length_y)<<"\t"<<
                          posZ<<"\t"<<
                          //fPlx_x<<"\t"<<fPlx_y<<"\t"<<fEdep_plx<<"\t"<<

                          endl;
                //  file_edep_zero_all_electron_arrived_detector<<edep/keV<<"\t"<<ekin/MeV<<"\t"<<procName<<"\t"<<particleName<<"\t"<<ekinpost<<"\t"<<ekinpost<<"\t"<<feventID<<endl;
        //file_edep_zero_all_electron_arrived_detector<<edep/keV<<"\t"<<ekin/MeV<<"\t"<<fEventAction->fEdep_plx<<"\t"<<fEventAction->fPlx_x<<"\t"<<fEventAction->fPlx_y<<"\t"<<track_id<<"\t"<<feventID<<endl;
            }
            file_edep_zero_all_electron_arrived_detector.close();
        }
        ///////////////////////////////////////////////
        if(particleName == "anti_nu_e"  && edep==0.) {
        fEventAction->AddE_zero_all_anti_nu_e_arrived_detector(edep, ekinpost, stepLength);///no dexdx
        G4String  edep_zero_all_anti_nu_e_arrived_detector;
        edep_zero_all_anti_nu_e_arrived_detector += fileName;
        edep_zero_all_anti_nu_e_arrived_detector +="_edep_zero_all_anti_nu_e_arrived_detector";
        edep_zero_all_anti_nu_e_arrived_detector +=".csv";
        ofstream file_edep_zero_all_anti_nu_e_arrived_detector (edep_zero_all_anti_nu_e_arrived_detector, ios::out | ios::app);
        if (file_edep_zero_all_anti_nu_e_arrived_detector.is_open()) {
        file_edep_zero_all_anti_nu_e_arrived_detector<<
                  trunc((pix_x))<<"\t"<<
                  trunc((pix_y))<<"\t"<<
                  edep/keV<<"\t"<<
                  ekinpost/MeV<<"\t"<<
                  ekinpre/MeV<<"\t"<<
                  edep/ekinpre*100<<"\t"<<
                  dEdL<<"\t"<<
                  feventID<<"\t"<<
                  //posX<<"\t"<<posY<<"\t"<<posZ<<"\t"<<
                  //round((pix_x))<<"\t"<<
                  //round((pix_y))<<"\t"<<
                  procName<<"\t"<<
                  parent_id<<"\t"<<
                  track_id<<"\t"<<

                  (posX+0.5*length_x)<<"\t"<<
                  (posY+0.5*length_y)<<"\t"<<
                  posZ<<"\t"<<
                  //fPlx_x<<"\t"<<fPlx_y<<"\t"<<fEdep_plx<<"\t"<<

                  endl;
        //  file_edep_zero_all_anti_nu_e_arrived_detector<<edep/keV<<"\t"<<ekin/MeV<<"\t"<<procName<<"\t"<<particleName<<"\t"<<ekinpost<<"\t"<<ekinpost<<"\t"<<feventID<<endl;
        //file_edep_zero_all_anti_nu_e_arrived_detector<<edep/keV<<"\t"<<ekin/MeV<<"\t"<<fEventAction->fEdep_plx<<"\t"<<fEventAction->fPlx_x<<"\t"<<fEventAction->fPlx_y<<"\t"<<track_id<<"\t"<<feventID<<endl;
        }
        file_edep_zero_all_anti_nu_e_arrived_detector.close();
        }
        ///////////////////////////////////////////////
        if(particleName == "alpha"  && edep==0.) {
            fEventAction->AddE_zero_all_alpha_arrived_detector(edep, ekinpost, stepLength);///no dexdx
            G4String  edep_zero_all_alpha_arrived_detector;
            edep_zero_all_alpha_arrived_detector += fileName;
            edep_zero_all_alpha_arrived_detector +="_edep_zero_all_alpha_arrived_detector";
            edep_zero_all_alpha_arrived_detector +=".csv";
            ofstream file_edep_zero_all_alpha_arrived_detector (edep_zero_all_alpha_arrived_detector, ios::out | ios::app);
            if (file_edep_zero_all_alpha_arrived_detector.is_open()) {
                file_edep_zero_all_alpha_arrived_detector<<
                          trunc((pix_x))<<"\t"<<
                          trunc((pix_y))<<"\t"<<
                          edep/keV<<"\t"<<
                          ekinpost/MeV<<"\t"<<
                          ekinpre/MeV<<"\t"<<
                          edep/ekinpre*100<<"\t"<<
                          dEdL<<"\t"<<
                          feventID<<"\t"<<
                          //posX<<"\t"<<posY<<"\t"<<posZ<<"\t"<<
                          //round((pix_x))<<"\t"<<
                          //round((pix_y))<<"\t"<<
                          procName<<"\t"<<
                          parent_id<<"\t"<<
                          track_id<<"\t"<<

                          (posX+0.5*length_x)<<"\t"<<
                          (posY+0.5*length_y)<<"\t"<<
                          posZ<<"\t"<<
                          //fPlx_x<<"\t"<<fPlx_y<<"\t"<<fEdep_plx<<"\t"<<

                          endl;
                //  file_edep_zero_all_alpha_arrived_detector<<edep/keV<<"\t"<<ekin/MeV<<"\t"<<procName<<"\t"<<particleName<<"\t"<<ekinpost<<"\t"<<ekinpost<<"\t"<<feventID<<endl;
        //file_edep_zero_all_alpha_arrived_detector<<edep/keV<<"\t"<<ekin/MeV<<"\t"<<fEventAction->fEdep_plx<<"\t"<<fEventAction->fPlx_x<<"\t"<<fEventAction->fPlx_y<<"\t"<<track_id<<"\t"<<feventID<<endl;
            }
            file_edep_zero_all_alpha_arrived_detector.close();
        }
        ///////////////////////////////////////////////
        if(particleName == "gamma"  && edep==0.) {
            fEventAction->AddE_zero_all_gamma_arrived_detector(edep, ekinpost, stepLength);///no dexdx
            G4String  edep_zero_all_gamma_arrived_detector;
            edep_zero_all_gamma_arrived_detector += fileName;
            edep_zero_all_gamma_arrived_detector +="_edep_zero_all_gamma_arrived_detector";
            edep_zero_all_gamma_arrived_detector +=".csv";
            ofstream file_edep_zero_all_gamma_arrived_detector (edep_zero_all_gamma_arrived_detector, ios::out | ios::app);
            if (file_edep_zero_all_gamma_arrived_detector.is_open()) {
                        file_edep_zero_all_gamma_arrived_detector<<
                        trunc((pix_x))<<"\t"<<
                        trunc((pix_y))<<"\t"<<
                        edep/keV<<"\t"<<
                        ekinpost/MeV<<"\t"<<
                        ekinpre/MeV<<"\t"<<
                        edep/ekinpre*100<<"\t"<<
                        dEdL<<"\t"<<
                        feventID<<"\t"<<
                        //posX<<"\t"<<posY<<"\t"<<posZ<<"\t"<<
                        //round((pix_x))<<"\t"<<
                        //round((pix_y))<<"\t"<<
                        procName<<"\t"<<
                        parent_id<<"\t"<<
                        track_id<<"\t"<<

                        (posX+0.5*length_x)<<"\t"<<
                        (posY+0.5*length_y)<<"\t"<<
                        posZ<<"\t"<<
                        //fPlx_x<<"\t"<<fPlx_y<<"\t"<<fEdep_plx<<"\t"<<

                        endl;
                //  file_edep_zero_all_gamma_arrived_detector<<edep/keV<<"\t"<<ekin/MeV<<"\t"<<procName<<"\t"<<particleName<<"\t"<<ekinpost<<"\t"<<ekinpost<<"\t"<<feventID<<endl;
        //file_edep_zero_all_gamma_arrived_detector<<edep/keV<<"\t"<<ekin/MeV<<"\t"<<fEventAction->fEdep_plx<<"\t"<<fEventAction->fPlx_x<<"\t"<<fEventAction->fPlx_y<<"\t"<<track_id<<"\t"<<feventID<<endl;
            }
            file_edep_zero_all_gamma_arrived_detector.close();
        }
        ////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////
        if(edep>0. ) {
            //fEventAction->AddE_dif0_all_particle_arrived_detector(edep, ekinpost);
            G4String  edep_dif0_all_particle_arrived_detector;
            edep_dif0_all_particle_arrived_detector += fileName;
            edep_dif0_all_particle_arrived_detector +="_edep_dif0_all_particle_arrived_detector";
            edep_dif0_all_particle_arrived_detector +=".csv";
            ofstream file_edep_dif0_all_particle_arrived_detector (edep_dif0_all_particle_arrived_detector, ios::out | ios::app);
            if (file_edep_dif0_all_particle_arrived_detector.is_open()) {
                file_edep_dif0_all_particle_arrived_detector<<
                          trunc((pix_x))<<"\t"<<
                          trunc((pix_y))<<"\t"<<
                          edep/keV<<"\t"<<
                          ekinpost/MeV<<"\t"<<
                          ekinpre<<"\t"<<
                          ekin<<"\t"<<
                          dEdL<<"\t"<<
                          feventID<<"\t"<<
                          particleName<<"\t"<<
                          procName<<"\t"<<
                          parent_id<<"\t"<<
                          track_id<<"\t"<<
                          posX<<"\t"<<
                          posY<<"\t"<<
                          posZ<<"\t"<<
                          endl;

            }
            file_edep_dif0_all_particle_arrived_detector.close();
        }
        //////////////////////////////////////////////
        if(particleName == "e-"  && edep>0.) {
            fEventAction->AddE_dif0_all_electron_arrived_detector(edep, ekinpost, stepLength);///no dexdx
            G4String  edep_dif0_all_electron_arrived_detector;
            edep_dif0_all_electron_arrived_detector += fileName;
            edep_dif0_all_electron_arrived_detector +="_edep_dif0_all_electron_arrived_detector";
            edep_dif0_all_electron_arrived_detector +=".csv";
            ofstream file_edep_dif0_all_electron_arrived_detector (edep_dif0_all_electron_arrived_detector, ios::out | ios::app);
            if (file_edep_dif0_all_electron_arrived_detector.is_open()) {
                file_edep_dif0_all_electron_arrived_detector<<
                          trunc((pix_x))<<"\t"<<
                          trunc((pix_y))<<"\t"<<
                          edep/keV<<"\t"<<
                          ekinpost/MeV<<"\t"<<
                          ekinpre<<"\t"<<
                          ekin<<"\t"<<
                          dEdL<<"\t"<<
                          feventID<<"\t"<<
                          //posX<<"\t"<<posY<<"\t"<<posZ<<"\t"<<
                          //round((pix_x))<<"\t"<<
                          //round((pix_y))<<"\t"<<
                          procName<<"\t"<<
                          parent_id<<"\t"<<
                          track_id<<"\t"<<

                          (posX+0.5*length_x)<<"\t"<<
                          (posY+0.5*length_y)<<"\t"<<
                          posZ<<"\t"<<
                          //fPlx_x<<"\t"<<fPlx_y<<"\t"<<fEdep_plx<<"\t"<<

                          endl;
                //  file_edep_dif0_all_electron_arrived_detector<<edep/keV<<"\t"<<ekin/MeV<<"\t"<<procName<<"\t"<<particleName<<"\t"<<ekinpost<<"\t"<<ekinpost<<"\t"<<feventID<<endl;
        //file_edep_dif0_all_electron_arrived_detector<<edep/keV<<"\t"<<ekin/MeV<<"\t"<<fEventAction->fEdep_plx<<"\t"<<fEventAction->fPlx_x<<"\t"<<fEventAction->fPlx_y<<"\t"<<track_id<<"\t"<<feventID<<endl;
            }
            file_edep_dif0_all_electron_arrived_detector.close();
        }
        ///////////////////////////////////////////////
        if(particleName == "anti_nu_e"  && edep>0.) {
        fEventAction->AddE_dif0_all_anti_nu_e_arrived_detector(edep, ekinpost, stepLength);///no dexdx
        G4String  edep_dif0_all_anti_nu_e_arrived_detector;
        edep_dif0_all_anti_nu_e_arrived_detector += fileName;
        edep_dif0_all_anti_nu_e_arrived_detector +="_edep_dif0_all_anti_nu_e_arrived_detector";
        edep_dif0_all_anti_nu_e_arrived_detector +=".csv";
        ofstream file_edep_dif0_all_anti_nu_e_arrived_detector (edep_dif0_all_anti_nu_e_arrived_detector, ios::out | ios::app);
        if (file_edep_dif0_all_anti_nu_e_arrived_detector.is_open()) {
        file_edep_dif0_all_anti_nu_e_arrived_detector<<
                  trunc((pix_x))<<"\t"<<
                  trunc((pix_y))<<"\t"<<
                  edep/keV<<"\t"<<
                  ekinpost/MeV<<"\t"<<
                  ekinpre/MeV<<"\t"<<
                  edep/ekinpre*100<<"\t"<<
                  dEdL<<"\t"<<
                  feventID<<"\t"<<
                  //posX<<"\t"<<posY<<"\t"<<posZ<<"\t"<<
                  //round((pix_x))<<"\t"<<
                  //round((pix_y))<<"\t"<<
                  procName<<"\t"<<
                  parent_id<<"\t"<<
                  track_id<<"\t"<<

                  (posX+0.5*length_x)<<"\t"<<
                  (posY+0.5*length_y)<<"\t"<<
                  posZ<<"\t"<<
                  //fPlx_x<<"\t"<<fPlx_y<<"\t"<<fEdep_plx<<"\t"<<

                  endl;
        //  file_edep_dif0_all_anti_nu_e_arrived_detector<<edep/keV<<"\t"<<ekin/MeV<<"\t"<<procName<<"\t"<<particleName<<"\t"<<ekinpost<<"\t"<<ekinpost<<"\t"<<feventID<<endl;
        //file_edep_dif0_all_anti_nu_e_arrived_detector<<edep/keV<<"\t"<<ekin/MeV<<"\t"<<fEventAction->fEdep_plx<<"\t"<<fEventAction->fPlx_x<<"\t"<<fEventAction->fPlx_y<<"\t"<<track_id<<"\t"<<feventID<<endl;
        }
        file_edep_dif0_all_anti_nu_e_arrived_detector.close();
        }
        ///////////////////////////////////////////////
        if(particleName == "alpha"  && edep>0.) {
            fEventAction->AddE_dif0_all_alpha_arrived_detector(edep, ekinpost, stepLength);///no dexdx
            G4String  edep_dif0_all_alpha_arrived_detector;
            edep_dif0_all_alpha_arrived_detector += fileName;
            edep_dif0_all_alpha_arrived_detector +="_edep_dif0_all_alpha_arrived_detector";
            edep_dif0_all_alpha_arrived_detector +=".csv";
            ofstream file_edep_dif0_all_alpha_arrived_detector (edep_dif0_all_alpha_arrived_detector, ios::out | ios::app);
            if (file_edep_dif0_all_alpha_arrived_detector.is_open()) {
                file_edep_dif0_all_alpha_arrived_detector<<
                          trunc((pix_x))<<"\t"<<
                          trunc((pix_y))<<"\t"<<
                          edep/keV<<"\t"<<
                          ekinpost/MeV<<"\t"<<
                          ekinpre/MeV<<"\t"<<
                          edep/ekinpre*100<<"\t"<<
                          dEdL<<"\t"<<
                          feventID<<"\t"<<
                          //posX<<"\t"<<posY<<"\t"<<posZ<<"\t"<<
                          //round((pix_x))<<"\t"<<
                          //round((pix_y))<<"\t"<<
                          procName<<"\t"<<
                          parent_id<<"\t"<<
                          track_id<<"\t"<<

                          (posX+0.5*length_x)<<"\t"<<
                          (posY+0.5*length_y)<<"\t"<<
                          posZ<<"\t"<<
                          //fPlx_x<<"\t"<<fPlx_y<<"\t"<<fEdep_plx<<"\t"<<

                          endl;
                //  file_edep_dif0_all_alpha_arrived_detector<<edep/keV<<"\t"<<ekin/MeV<<"\t"<<procName<<"\t"<<particleName<<"\t"<<ekinpost<<"\t"<<ekinpost<<"\t"<<feventID<<endl;
        //file_edep_dif0_all_alpha_arrived_detector<<edep/keV<<"\t"<<ekin/MeV<<"\t"<<fEventAction->fEdep_plx<<"\t"<<fEventAction->fPlx_x<<"\t"<<fEventAction->fPlx_y<<"\t"<<track_id<<"\t"<<feventID<<endl;
            }
            file_edep_dif0_all_alpha_arrived_detector.close();
        }
        ///////////////////////////////////////////////
        if(particleName == "gamma"  && edep>0.) {
            fEventAction->AddE_dif0_all_gamma_arrived_detector(edep, ekinpost, stepLength);///no dexdx
            G4String  edep_dif0_all_gamma_arrived_detector;
            edep_dif0_all_gamma_arrived_detector += fileName;
            edep_dif0_all_gamma_arrived_detector +="_edep_dif0_all_gamma_arrived_detector";
            edep_dif0_all_gamma_arrived_detector +=".csv";
            ofstream file_edep_dif0_all_gamma_arrived_detector (edep_dif0_all_gamma_arrived_detector, ios::out | ios::app);
            if (file_edep_dif0_all_gamma_arrived_detector.is_open()) {
                        file_edep_dif0_all_gamma_arrived_detector<<
                        trunc((pix_x))<<"\t"<<
                        trunc((pix_y))<<"\t"<<
                        edep/keV<<"\t"<<
                        ekinpost/MeV<<"\t"<<
                        ekinpre/MeV<<"\t"<<
                        edep/ekinpre*100<<"\t"<<
                        dEdL<<"\t"<<
                        feventID<<"\t"<<
                        //posX<<"\t"<<posY<<"\t"<<posZ<<"\t"<<
                        //round((pix_x))<<"\t"<<
                        //round((pix_y))<<"\t"<<
                        procName<<"\t"<<
                        parent_id<<"\t"<<
                        track_id<<"\t"<<

                        (posX+0.5*length_x)<<"\t"<<
                        (posY+0.5*length_y)<<"\t"<<
                        posZ<<"\t"<<
                        //fPlx_x<<"\t"<<fPlx_y<<"\t"<<fEdep_plx<<"\t"<<

                        endl;
                //  file_edep_dif0_all_gamma_arrived_detector<<edep/keV<<"\t"<<ekin/MeV<<"\t"<<procName<<"\t"<<particleName<<"\t"<<ekinpost<<"\t"<<ekinpost<<"\t"<<feventID<<endl;
        //file_edep_dif0_all_gamma_arrived_detector<<edep/keV<<"\t"<<ekin/MeV<<"\t"<<fEventAction->fEdep_plx<<"\t"<<fEventAction->fPlx_x<<"\t"<<fEventAction->fPlx_y<<"\t"<<track_id<<"\t"<<feventID<<endl;
            }
            file_edep_dif0_all_gamma_arrived_detector.close();
        }
        ////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////

        /// PRIMARIES //
        if (track_id ==1 && parent_id == 0) {
            //  std::cout<<"Primary_Particle_in_Detectror: "<<particleName<<" Edep_step(MeV): "<<edep/keV<<std::endl;
            //std::cout<<"Track #"<<track_id<<" of "<<particleName<<" E_step(MeV)= "<<edep/keV<<" produced by Track ID= "<<parent_id<<std::endl;
            if(particleName == "e-" ) {
                fEventAction->AddE_electron_primary_detector(edep, ekinpost);
                //fEventAction->AddKe_electron_primary_detector(edep);
                G4String  edep_electron_primary_detector;
                edep_electron_primary_detector += fileName;
                edep_electron_primary_detector +="_edep_electron_primary_detector";
                edep_electron_primary_detector +=".csv";
                ofstream file_edep_electron_primary_detector (edep_electron_primary_detector, ios::out | ios::app);
                if (file_edep_electron_primary_detector.is_open()) {
                    file_edep_electron_primary_detector<<
                              edep/keV<<"\t"<<
                              ekinpost/MeV<<"\t"<<
                              //procName<<"\t"<<
                              parent_id<<"\t"<<
                              track_id<<"\t"<<
feventID<<endl;
                }
                file_edep_electron_primary_detector.close();
            }
        }
        /// SECONDARIES ///
        if ( track_id != 1 && parent_id > 0 ) {
            //  std::cout<<"Secondary_Particle_in_Detectror: "<<particleName<<" Edep_step(MeV): "<<edep/keV
            //         <<" Track_# "<<track_id<<" produced_by_Track_ID= "<<parent_id
            //       <<std::endl;
            //std::cout<<"Track #"<<track_id<<" of "<<particleName<<" E_step(MeV)= "<<edep/keV<<" produced by Track ID= "<<parent_id<<std::endl;
            if(particleName == "e-" ) {
                fEventAction->AddE_electron_second_detector(edep, ekinpost);
                G4String  edep_electron_second_detector;
                edep_electron_second_detector += fileName;
                edep_electron_second_detector +="_edep_electron_second_detector";
                edep_electron_second_detector +=".csv";
                ofstream file_edep_electron_second_detector (edep_electron_second_detector, ios::out | ios::app);
                if (file_edep_electron_second_detector.is_open()) {
                    file_edep_electron_second_detector<<
                              edep/keV<<"\t"<<
                              ekinpost/MeV<<"\t"<<
                              //procName<<"\t"<<
                              parent_id<<"\t"<<
                              track_id<<"\t"<<
feventID<<endl;
                }
                file_edep_electron_second_detector.close();

            }
            if(particleName == "gamma" ) {
                fEventAction->AddE_gamma_second_detector(edep, ekinpost);
                G4String  edep_gamma_second_detector;
                edep_gamma_second_detector += fileName;
                edep_gamma_second_detector +="_edep_gamma_second_detector";
                edep_gamma_second_detector +=".csv";
                ofstream file_edep_gamma_second_detector (edep_gamma_second_detector, ios::out | ios::app);
                if (file_edep_gamma_second_detector.is_open()) {
                    file_edep_gamma_second_detector<<
                              edep/keV<<"\t"<<
                              ekinpost/MeV<<"\t"<<
                              //procName<<"\t"<<
                              parent_id<<"\t"<<
                              track_id<<"\t"<<
feventID<<endl;
                }
                file_edep_gamma_second_detector.close();

            }
        }


        ////////electron Particle_is_a_Secondary_but_Parent_was_elec
        if ( (track_id != 1) && (parent_id == 4) ) {
            //std::cout<<"Particle_is_a_Secondary_but_Parent_was_elec_in_Detectror: "<<particleName<<" Edep_step(MeV): "<<edep/keV<<std::endl;
            //std::cout<<"Track #"<<track_id<<" of "<<particleName<<" E_step(MeV)= "<<edep/keV<<" produced by Track ID= "<<parent_id<<std::endl;
            if(particleName == "e-" ) {
                fEventAction->AddE_electron_second_parent_elec_detector(edep, ekinpost);
                G4String  edep_electron_second_parent_elec_detector;
                edep_electron_second_parent_elec_detector += fileName;
                edep_electron_second_parent_elec_detector +="_edep_electron_second_parent_elec_detector";
                edep_electron_second_parent_elec_detector +=".csv";
                ofstream file_edep_electron_second_parent_elec_detector (edep_electron_second_parent_elec_detector, ios::out | ios::app);
                if (file_edep_electron_second_parent_elec_detector.is_open()) {
                    file_edep_electron_second_parent_elec_detector<<
                              edep/keV<<"\t"<<
                              ekinpost/MeV<<"\t"<<
                              //procName<<"\t"<<
                              parent_id<<"\t"<<
                              track_id<<"\t"<<
feventID<<endl;
                }
                file_edep_electron_second_parent_elec_detector.close();

            }
        }



        if(fSource =="e" || fSource =="alpha") {
        }
  else{
          //G4cout << "testttttttttttttttttttttttttttt: "  << G4endl;


        ////////electron_second_parent_Sr90_detector
        if( track_id == 4 && parent_id == 1 ) {
            ////////electron_second_parent_Sr90_detector
            if(particleName == "e-" ) {
                fEventAction->AddE_electron_second_parent_Sr90_detector(edep, ekinpost);
                //G4cout<<" Process_Name: "<<procName<<" Particles_Name: "<<particleName<<G4endl;
                G4String  edep_electron_second_parent_Sr90_detector;
                edep_electron_second_parent_Sr90_detector += fileName;
                edep_electron_second_parent_Sr90_detector +="_edep_electron_second_parent_Sr90_detector";
                edep_electron_second_parent_Sr90_detector +=".csv";
                ofstream file_edep_electron_second_parent_Sr90_detector ( edep_electron_second_parent_Sr90_detector, ios::out | ios::app);
                if (file_edep_electron_second_parent_Sr90_detector.is_open()) {
                    file_edep_electron_second_parent_Sr90_detector<<
                              edep/keV<<"\t"<<
                              ekinpost/MeV<<"\t"<<
feventID<<endl;
                              //"\t"<<procName<<endl;
                }
                file_edep_electron_second_parent_Sr90_detector.close();
            }
            ////////electron_second_transport_parent_Sr90_detector
            if(particleName == "e-" && procName == "Transportation" ) {
                fEventAction->AddE_electron_second_transport_parent_Sr90_detector(edep, ekinpost);
                //G4cout<<" Process_Name: "<<procName<<" Particles_Name: "<<particleName<<G4endl;
                G4String  edep_electron_second_transport_parent_Sr90_detector;
                edep_electron_second_transport_parent_Sr90_detector += fileName;
                edep_electron_second_transport_parent_Sr90_detector +="_edep_electron_second_transport_parent_Sr90_detector";
                edep_electron_second_transport_parent_Sr90_detector +=".csv";
                ofstream file_edep_electron_second_transport_parent_Sr90_detector ( edep_electron_second_transport_parent_Sr90_detector, ios::out | ios::app);
                if (file_edep_electron_second_transport_parent_Sr90_detector.is_open()) {
                    file_edep_electron_second_transport_parent_Sr90_detector<<
                            edep/keV<<"\t"<<ekinpost/MeV<<"\t"<<
feventID<<endl;
                }
                file_edep_electron_second_transport_parent_Sr90_detector.close();
            }
            ////////electron_second_eIoni_parent_Sr90_detector
            if(particleName == "e-" && procName == "eIoni" ) {
                fEventAction->AddE_electron_second_eIoni_parent_Sr90_detector(edep, ekinpost);
                //G4cout<<" Process_Name: "<<procName<<" Particles_Name: "<<particleName<<G4endl;
                G4String  edep_electron_second_eIoni_parent_Sr90_detector;
                edep_electron_second_eIoni_parent_Sr90_detector += fileName;
                edep_electron_second_eIoni_parent_Sr90_detector +="_edep_electron_second_eIoni_parent_Sr90_detector";
                edep_electron_second_eIoni_parent_Sr90_detector +=".csv";
                ofstream file_edep_electron_second_eIoni_parent_Sr90_detector ( edep_electron_second_eIoni_parent_Sr90_detector, ios::out | ios::app);
                if (file_edep_electron_second_eIoni_parent_Sr90_detector.is_open()) {
                    file_edep_electron_second_eIoni_parent_Sr90_detector<<
                            edep/keV<<"\t"<<ekinpost/MeV<<"\t"<<
feventID<<endl;
                }
                file_edep_electron_second_eIoni_parent_Sr90_detector.close();
            }
            ////////electron_second_eBrem_parent_Sr90_detector
            if(particleName == "e-" && procName == "eBrema" ) {
                fEventAction->AddE_electron_second_eBrem_parent_Sr90_detector(edep, ekinpost);
                //G4cout<<" Process_Name: "<<procName<<" Particles_Name: "<<particleName<<G4endl;
                G4String  edep_electron_second_eBrem_parent_Sr90_detector;
                edep_electron_second_eBrem_parent_Sr90_detector += fileName;
                edep_electron_second_eBrem_parent_Sr90_detector +="_edep_electron_second_eBrem_parent_Sr90_detector";
                edep_electron_second_eBrem_parent_Sr90_detector +=".csv";
                ofstream file_edep_electron_second_eBrem_parent_Sr90_detector ( edep_electron_second_eBrem_parent_Sr90_detector, ios::out | ios::app);
                if (file_edep_electron_second_eBrem_parent_Sr90_detector.is_open()) {
                    file_edep_electron_second_eBrem_parent_Sr90_detector<<
                            edep/keV<<"\t"<<ekinpost/MeV<<"\t"<<
feventID<<endl;
                }
                file_edep_electron_second_eBrem_parent_Sr90_detector.close();
            }
            ////////electron_second_msc_parent_Sr90_detector
            if(particleName == "e-" && procName == "msc" ) {
                fEventAction->AddE_electron_second_msc_parent_Sr90_detector(edep, ekinpost);
                //G4cout<<" Process_Name: "<<procName<<" Particles_Name: "<<particleName<<G4endl;
                G4String  edep_electron_second_msc_parent_Sr90_detector;
                edep_electron_second_msc_parent_Sr90_detector += fileName;
                edep_electron_second_msc_parent_Sr90_detector +="_edep_electron_second_msc_parent_Sr90_detector";
                edep_electron_second_msc_parent_Sr90_detector +=".csv";
                ofstream file_edep_electron_second_msc_parent_Sr90_detector ( edep_electron_second_msc_parent_Sr90_detector, ios::out | ios::app);
                if (file_edep_electron_second_msc_parent_Sr90_detector.is_open()) {
                    file_edep_electron_second_msc_parent_Sr90_detector<<
                            edep/keV<<"\t"<<ekinpost/MeV<<"\t"<<
feventID<<endl;
                }
                file_edep_electron_second_msc_parent_Sr90_detector.close();
            }

        }


        ////////Y90 Particle_is_a_Secondary_but_Parent_was_Y90
        if ( (track_id != 1) && (parent_id == 2) ) {
            //std::cout<<"Particle_is_a_Secondary_but_Parent_was_Y90_in_Detectror: "<<particleName<<" Edep_step(MeV): "<<edep/keV<<std::endl;
            //std::cout<<"Track #"<<track_id<<" of "<<particleName<<" E_step(MeV)= "<<edep/keV<<" produced by Track ID= "<<parent_id<<std::endl;
            if(particleName == "e-" ) {
                fEventAction->AddE_electron_second_parent_Y90_detector(edep, ekinpost);
                G4String  edep_electron_second_parent_Y90_detector;
                edep_electron_second_parent_Y90_detector += fileName;
                edep_electron_second_parent_Y90_detector +="_edep_electron_second_parent_Y90_detector";
                edep_electron_second_parent_Y90_detector +=".csv";
                ofstream file_edep_electron_second_parent_Y90_detector (edep_electron_second_parent_Y90_detector, ios::out | ios::app);
                if (file_edep_electron_second_parent_Y90_detector.is_open()) {
                    file_edep_electron_second_parent_Y90_detector<<
                              edep/keV<<"\t"<<
                              ekinpost/MeV<<"\t"<<
                              //procName<<"\t"<<
                              parent_id<<"\t"<<
                              track_id<<"\t"<<
feventID<<endl;
                }
                file_edep_electron_second_parent_Y90_detector.close();

            }
        }


/////////////////////////////////////////////////
      }
///////////////////////////////////////////////////
    }

////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
