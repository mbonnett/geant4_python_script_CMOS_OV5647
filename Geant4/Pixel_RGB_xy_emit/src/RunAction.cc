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
/// \file RunAction.cc
/// \brief Implementation of the RunAction class

#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "RunActionMessenger.hh"

//#include "g4root.hh"
///#include "Analysis.hh"
//#include "G4AnalysisManager.hh"

#include "Randomize.hh"
#include <ctime>

#include <iostream>
#include <fstream>
using namespace std;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
    : G4UserRunAction(),
      fRunActionMessenger(0),
//fFilename("Pixel_Si"),
      fGoodEvents(0)

{
    fSource = "gamma";
    fThreshold = 0*eV;
    fFilename = "Pixel_Si";

    fN_pixel_x = 2592;
    fN_pixel_y = 1944;
    fpixel_x = 1.4*um;
    fpixel_y = 1.4*um;

    fSaveProcPart = true;
    // Register accumulable to the accumulable manager
    G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
    accumulableManager->RegisterAccumulable(fGoodEvents);

    fRunActionMessenger = new RunActionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
    //  delete G4AnalysisManager::Instance();
    delete fRunActionMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* run)
{
    G4cout<< "### Run "<< run->GetRunID()<< " start."<< G4endl;

    G4String  fileName;
    fileName = fFilename;
    fileName += "_";
    fileName += fSource;

    G4String pxledep2Dname;
    pxledep2Dname += fileName;
    pxledep2Dname +="_pixel_edep2D";
    pxledep2Dname +=".csv";
    ofstream file_pxledep2D (pxledep2Dname, ios::out | ios::trunc);
    file_pxledep2D<<"x("<<fN_pixel_x<<")\t"<<"y("<<fN_pixel_y<<")\t"<<"E_dep(keV)"<<"\t"//<<"E_emit(keV)"<<"\t"
    <<"eventID"<<endl;


    G4String pxl_Edep0emit_name;
    pxl_Edep0emit_name += fileName;
    pxl_Edep0emit_name +="_pixel_edep_zero_emit.csv";
    ofstream pxl_Edep0emit (pxl_Edep0emit_name, ios::out | ios::trunc);
    pxl_Edep0emit<<"E_dep(keV)"<<"\t"<<"E_emit(keV)"<<"\t"<<"eventID"<<endl;


    G4String pxl_Edepemit_name;
    pxl_Edepemit_name += fileName;
    pxl_Edepemit_name +="_pixel_edep_emit.csv";
    ofstream pxl_Edepemit (pxl_Edepemit_name, ios::out | ios::trunc);
    pxl_Edepemit<<"E_dep(keV)"<<"\t"<<"E_emit(keV)"<<"\t"<<"eventID"<<endl;



//   file_pxledep<<"x"<<"\t"<<"y"<<"\t"<<"Edep(keV)"<<endl;//
  //  file_pxledep2D<<fN_pixel_x<<"\t"<<fN_pixel_y<<"\t"<<0<<endl;

/*
    G4String  edep_pxlname;
    edep_pxlname += fileName;
    edep_pxlname +="_edep_pixel";
    edep_pxlname +=".csv";
    ofstream file_edep_pxl (edep_pxlname, ios::out | ios::trunc);
    file_edep_pxl<<Edep(keV)"<<"\t"<<
"eventID"<<endl;
*/
    G4String  edep_detcname;
    edep_detcname += fileName;
    edep_detcname += "_edep_dect";
    edep_detcname += ".csv";
    ofstream file_edep_detc (edep_detcname, ios::out | ios::trunc);
    file_edep_detc<<"E(keV)"<<"\t"<<
  "eventID"<<endl;


/*
    G4String  hit_pxl2Dname;
    hit_pxl2Dname += fileName;
    hit_pxl2Dname +="_hit_pixel2D";
    hit_pxl2Dname +=".csv";
    ofstream file_hit_pxl2D (hit_pxl2Dname, ios::out | ios::trunc);
    file_hit_pxl2D<<"x("<<fN_pixel_x<<")\t"<<"y("<<fN_pixel_y<<")\t"<<"hit"<<endl;
//   file_pxledep<<"x"<<"\t"<<"y"<<"\t"<<"Edep(keV)"<<endl;//
//    file_hit_pxl2D<<fN_pixel_x<<"\t"<<fN_pixel_y<<"\t"<<0<<endl;

    G4String  hit_pxlname;
    hit_pxlname += fileName;
    hit_pxlname +="_hit_pixel";
    hit_pxlname +=".csv";
    ofstream file_hit_pxl (hit_pxlname, ios::out | ios::trunc);
    file_hit_pxl<<"hit"<<endl;
*/

////////////////////////////////////////////////////////////////////////////////
if(fSource=="e" || fSource=="alpha") {
  G4String  emit_ParticlesCharge;
  emit_ParticlesCharge += fileName;
  emit_ParticlesCharge +="_emit_ParticlesCharge";
  emit_ParticlesCharge +=".csv";
  ofstream file_emit_ParticlesCharge (emit_ParticlesCharge, ios::out | ios::trunc);
  file_emit_ParticlesCharge<<"E(MeV)"<<"\t"<<"totalE"<<"\t"<<
"eventID"<<endl;

  G4String  emit_electron;
  emit_electron += fileName;
  emit_electron +="_emit_electron";
  emit_electron +=".csv";
  ofstream file_emit_electron (emit_electron, ios::out | ios::trunc);
  file_emit_electron<<"E(MeV)"<<"\t"<<"totalE"<<"\t"<<
"eventID"<<endl;

  G4String  emit_positron;
  emit_positron += fileName;
  emit_positron +="_emit_positron";
  emit_positron +=".csv";
  ofstream file_emit_positron (emit_positron, ios::out | ios::trunc);
  file_emit_positron<<"E(MeV)"<<"\t"<<"totalE"<<"\t"<<
"eventID"<<endl;

  G4String  emit_alpha;
  emit_alpha += fileName;
  emit_alpha +="_emit_alpha";
  emit_alpha +=".csv";
  ofstream file_emit_alpha (emit_alpha, ios::out | ios::trunc);
  file_emit_alpha<<"E(MeV)"<<"\t"<<"totalE"<<"\t"<<
"eventID"<<endl;

  G4String  emit_gamma;
  emit_gamma += fileName;
  emit_gamma +="_emit_gamma";
  emit_gamma +=".csv";
  ofstream file_emit_gamma (emit_gamma, ios::out | ios::trunc);
  file_emit_gamma<<"E(MeV)"<<"\t"<<"totalE"<<"\t"<<
"eventID"<<endl;

}
else{
    G4String  emit_allRDecayProducts;
    emit_allRDecayProducts += fileName;
    emit_allRDecayProducts +="_emit_allRDecayProducts";
    emit_allRDecayProducts +=".csv";
    ofstream file_emit_allRDecayProducts (emit_allRDecayProducts, ios::out | ios::trunc);
    file_emit_allRDecayProducts<<"E(MeV)"<<"\t"<<"totalE"<<"\t"<<"eventID"<<endl;;

    G4String  emit_ParticlesCharge;
    emit_ParticlesCharge += fileName;
    emit_ParticlesCharge +="_emit_ParticlesCharge";
    emit_ParticlesCharge +=".csv";
    ofstream file_emit_ParticlesCharge (emit_ParticlesCharge, ios::out | ios::trunc);
    file_emit_ParticlesCharge<<"E(MeV)"<<"\t"<<"totalE"<<"\t"<<"eventID"<<endl;;

    G4String  emit_electron;
    emit_electron += fileName;
    emit_electron +="_emit_electron";
    emit_electron +=".csv";
    ofstream file_emit_electron (emit_electron, ios::out | ios::trunc);
    file_emit_electron<<"E(MeV)"<<"\t"<<"totalE"<<"\t"<<"eventID"<<endl;;

    G4String  emit_positron;
    emit_positron += fileName;
    emit_positron +="_emit_positron";
    emit_positron +=".csv";
    ofstream file_emit_positron (emit_positron, ios::out | ios::trunc);
    file_emit_positron<<"E(MeV)"<<"\t"<<"totalE"<<"\t"<<"eventID"<<endl;;

    G4String  emit_alpha;
    emit_alpha += fileName;
    emit_alpha +="_emit_alpha";
    emit_alpha +=".csv";
    ofstream file_emit_alpha (emit_alpha, ios::out | ios::trunc);
    file_emit_alpha<<"E(MeV)"<<"\t"<<"totalE"<<"\t"<<"eventID"<<endl;;

    G4String  emit_gamma;
    emit_gamma += fileName;
    emit_gamma +="_emit_gamma";
    emit_gamma +=".csv";
    ofstream file_emit_gamma (emit_gamma, ios::out | ios::trunc);
    file_emit_gamma<<"E(MeV)"<<"\t"<<"totalE"<<"\t"<<"eventID"<<endl;;

    G4String  emit_anti_nu_elec;
    emit_anti_nu_elec += fileName;
    emit_anti_nu_elec +="_emit_anti_nu_elec";
    emit_anti_nu_elec +=".csv";
    ofstream file_emit_anti_nu_elec (emit_anti_nu_elec, ios::out | ios::trunc);
    file_emit_anti_nu_elec<<"E(MeV)"<<"\t"<<"totalE"<<"\t"<<"eventID"<<endl;;

    if(fSource=="Sr90") {
        G4String  emit_Y90;
        emit_Y90 += fileName;
        emit_Y90 +="_emit_Y90";
        emit_Y90 +=".csv";
        ofstream file_emit_Y90 (emit_Y90, ios::out | ios::trunc);
        file_emit_Y90<<"E(MeV)"<<"\t"<<"totalE"<<"\t"<<"eventID"<<endl;;

        G4String  emit_Zr90;
        emit_Zr90 += fileName;
        emit_Zr90 +="_emit_Zr90";
        emit_Zr90 +=".csv";
        ofstream file_emit_Zr90 (emit_Zr90, ios::out | ios::trunc);
        file_emit_Zr90<<"E(MeV)"<<"\t"<<"totalE"<<"\t"<<"eventID"<<endl;;
    }


    if(fSource=="Cs137") {
        G4String  emit_Ba137_661_659;
        emit_Ba137_661_659 += fileName;
        emit_Ba137_661_659 +="_emit_Ba137_661_659";
        emit_Ba137_661_659 +=".csv";
        ofstream file_emit_Ba137_661_659 (emit_Ba137_661_659, ios::out | ios::trunc);
        file_emit_Ba137_661_659<<"E(MeV)"<<"\t"<<"totalE"<<"\t"<<"eventID"<<endl;;

        G4String  emit_Ba137;
        emit_Ba137 += fileName;
        emit_Ba137 +="_emit_Ba137";
        emit_Ba137 +=".csv";
        ofstream file_emit_Ba137 (emit_Ba137, ios::out | ios::trunc);
        file_emit_Ba137<<"E(MeV)"<<"\t"<<"totalE"<<"\t"<<"eventID"<<endl;;
    }

    if(fSource=="Fe55") {
        G4String  emit_Mn55;
        emit_Mn55 += fileName;
        emit_Mn55 +="_emit_Mn55";
        emit_Mn55 +=".csv";
        ofstream file_emit_Mn55 (emit_Mn55, ios::out | ios::trunc);
        file_emit_Mn55<<"E(MeV)"<<"\t"<<"totalE"<<"\t"<<"eventID"<<endl;;
    }


   }
////////////////////////////////////////////////////////////////////////////////

G4String  charge_dep2Dname;
charge_dep2Dname += fileName;
charge_dep2Dname +="_charge_dep2D";
charge_dep2Dname +=".csv";
ofstream file_charge_dep2D (charge_dep2Dname, ios::out | ios::trunc);
file_charge_dep2D<<"x("<<fN_pixel_x<<")\t"<<"y("<<fN_pixel_y<<")\t"<<"Q"<<"\t"<<
"eventID"<<endl;

//  file_charge_dep2D<<fN_pixel_x<<"\t"<<fN_pixel_y<<"\t"<<0<<endl;

G4String    gamma_dep2Dname;
gamma_dep2Dname += fileName;
gamma_dep2Dname +="_gamma_dep2D";
gamma_dep2Dname +=".csv";
ofstream file_gamma_dep2D (gamma_dep2Dname, ios::out | ios::trunc);
file_gamma_dep2D<<"x("<<fN_pixel_x<<")\t"<<"y("<<fN_pixel_y<<")\t"<<"gamma"<<"\t"<<
"eventID"<<endl;

//  file_gamma_dep2D<<fN_pixel_x<<"\t"<<fN_pixel_y<<"\t"<<0<<endl;

G4String  electron_dep2Dname;
electron_dep2Dname += fileName;
electron_dep2Dname +="_electron_dep2D";
electron_dep2Dname +=".csv";
ofstream file_electron_dep2D  (electron_dep2Dname, ios::out | ios::trunc);
file_electron_dep2D <<"x("<<fN_pixel_x<<")\t"<<"y("<<fN_pixel_y<<")\t"<<"e-"<<"\t"<<
"eventID"<<endl;

//  file_electron_dep2D <<fN_pixel_x<<"\t"<<fN_pixel_y<<"\t"<<0<<endl


    G4String  alpha_edep2Dname;
    alpha_edep2Dname += fileName;
    alpha_edep2Dname +="_alpha_edep2D";
    alpha_edep2Dname +=".csv";
    ofstream file_alpha_edep2D (alpha_edep2Dname, ios::out | ios::trunc);
    file_alpha_edep2D<<"x("<<fN_pixel_x<<")\t"<<"y("<<fN_pixel_y<<")\t"<<"Edep(keV)"<<"\t"<<
  "eventID"<<endl;

  //  file_alpha_edep2D<<fN_pixel_x<<"\t"<<fN_pixel_y<<"\t"<<0<<endl;

  G4String  edep_alpha_detcname;
  edep_alpha_detcname += fileName;
  edep_alpha_detcname += "_alpha_edep_dect";
  edep_alpha_detcname += ".csv";
  ofstream file_alpha_edep_detc (edep_alpha_detcname, ios::out | ios::trunc);
  file_alpha_edep_detc<<"Edep(keV)"<<"\t"<<"eventID"<<endl;


    G4String    gamma_edep2Dname;
    gamma_edep2Dname += fileName;
    gamma_edep2Dname +="_gamma_edep2D";
    gamma_edep2Dname +=".csv";
    ofstream file_gamma_edep2D (gamma_edep2Dname, ios::out | ios::trunc);
    file_gamma_edep2D<<"x("<<fN_pixel_x<<")\t"<<"y("<<fN_pixel_y<<")\t"<<"Edep(keV)"<<"\t"<<
  "eventID"<<endl;

  //  file_gamma_edep2D<<fN_pixel_x<<"\t"<<fN_pixel_y<<"\t"<<0<<endl;

    G4String  electron_edep2Dname;
    electron_edep2Dname += fileName;
    electron_edep2Dname +="_electron_edep2D";
    electron_edep2Dname +=".csv";
    ofstream file_electron_edep2D (electron_edep2Dname, ios::out | ios::trunc);
    file_electron_edep2D<<"x("<<fN_pixel_x<<")\t"<<"y("<<fN_pixel_y<<")\t"<<"Edep(keV)"<<"\t"<<
  "eventID"<<endl;

  //  file_electron_edep2D<<fN_pixel_x<<"\t"<<fN_pixel_y<<"\t"<<0<<endl;

  G4String  charge_dep_detcname;
  charge_dep_detcname += fileName;
  charge_dep_detcname += "_charge_dep_dect";
  charge_dep_detcname += ".csv";
  ofstream file_charge_dep_detc (charge_dep_detcname, ios::out | ios::trunc);
  file_charge_dep_detc<<"q"<<"\t"<<
"eventID"<<endl;


  G4String  gamma_dep_detcname;
  gamma_dep_detcname += fileName;
  gamma_dep_detcname += "_gamma_dep_dect";
  gamma_dep_detcname += ".csv";
  ofstream file_gamma_dep_detc (gamma_dep_detcname, ios::out | ios::trunc);
  file_gamma_dep_detc<<"gamma"<<"\t"<<
"eventID"<<endl;


  G4String  electron_dep_detcname;
  electron_dep_detcname += fileName;
  electron_dep_detcname += "_electron_dep_dect";
  electron_dep_detcname += ".csv";
  ofstream file_electron_dep_detc (electron_dep_detcname, ios::out | ios::trunc);
  file_electron_dep_detc<<"electron"<<"\t"<<
"eventID"<<endl;



  if(fSource =="e" || fSource =="alpha") {
  }
  else{
  //  G4cout << "testttttttttttttttttttttttttttt: "  << G4endl;


    G4String  positron_edep2Dname;
    positron_edep2Dname += fileName;
    positron_edep2Dname +="_positron_edep2D";
    positron_edep2Dname +=".csv";
    ofstream file_positron_edep2D (positron_edep2Dname, ios::out | ios::trunc);
    file_positron_edep2D<<"x("<<fN_pixel_x<<")\t"<<"y("<<fN_pixel_y<<")\t"<<"Edep(keV)"<<"\t"<<
  "eventID"<<endl;

  //  file_positron_edep2D<<fN_pixel_x<<"\t"<<fN_pixel_y<<"\t"<<0<<endl;

    G4String  muonMinus_edep2Dname;
    muonMinus_edep2Dname += fileName;
    muonMinus_edep2Dname +="_muonMinus_edep2D";
    muonMinus_edep2Dname +=".csv";
    ofstream file_muonMinus_edep2D (muonMinus_edep2Dname, ios::out | ios::trunc);
    file_muonMinus_edep2D<<"x("<<fN_pixel_x<<")\t"<<"y("<<fN_pixel_y<<")\t"<<"Edep(keV)"<<"\t"<<
  "eventID"<<endl;

  //  file_muonMinus_edep2D<<fN_pixel_x<<"\t"<<fN_pixel_y<<"\t"<<0<<endl;

    G4String  muonPlus_edep2Dname;
    muonPlus_edep2Dname += fileName;
    muonPlus_edep2Dname +="_muonPlus_edep2D";
    muonPlus_edep2Dname +=".csv";
    ofstream file_muonPlus_edep2D (muonPlus_edep2Dname, ios::out | ios::trunc);
    file_muonPlus_edep2D<<"x("<<fN_pixel_x<<")\t"<<"y("<<fN_pixel_y<<")\t"<<"Edep(keV)"<<"\t"<<
  "eventID"<<endl;

  //  file_muonPlus_edep2D<<fN_pixel_x<<"\t"<<fN_pixel_y<<"\t"<<0<<endl;

///////////////////////////////////////////////////////////////////////////////


    G4String  edep_gamma_detcname;
    edep_gamma_detcname += fileName;
    edep_gamma_detcname += "_gamma_edep_dect";
    edep_gamma_detcname += ".csv";
    ofstream file_gamma_edep_detc (edep_gamma_detcname, ios::out | ios::trunc);
    file_gamma_edep_detc<<"E(keV)"<<"\t"<<
  "eventID"<<endl;


    G4String  edep_electron_detcname;
    edep_electron_detcname += fileName;
    edep_electron_detcname += "_electron_edep_dect";
    edep_electron_detcname += ".csv";
    ofstream file_electron_edep_detc (edep_electron_detcname, ios::out | ios::trunc);
    file_electron_edep_detc<<"E(keV)"<<"\t"<<
  "eventID"<<endl;


    G4String  edep_positron_detcname;
    edep_positron_detcname += fileName;
    edep_positron_detcname += "_positron_edep_dect";
    edep_positron_detcname += ".csv";
    ofstream file_positron_edep_detc (edep_positron_detcname, ios::out | ios::trunc);
    file_positron_edep_detc<<"E(keV)"<<"\t"<<
  "eventID"<<endl;


    G4String  edep_alpha_detcname;
    edep_alpha_detcname += fileName;
    edep_alpha_detcname += "_alpha_edep_dect";
    edep_alpha_detcname += ".csv";
    ofstream file_alpha_edep_detc (edep_alpha_detcname, ios::out | ios::trunc);
    file_alpha_edep_detc<<"E(keV)"<<"\t"<<
  "eventID"<<endl;


    G4String  edep_muonMinus_detcname;
    edep_muonMinus_detcname += fileName;
    edep_muonMinus_detcname += "_muonMinus_edep_dect";
    edep_muonMinus_detcname += ".csv";
    ofstream file_muonMinus_edep_detc (edep_muonMinus_detcname, ios::out | ios::trunc);
    file_muonMinus_edep_detc<<"E(keV)"<<"\t"<<
  "eventID"<<endl;


    G4String  edep_muonPlus_detcname;
    edep_muonPlus_detcname += fileName;
    edep_muonPlus_detcname += "_muonPlus_edep_dect";
    edep_muonPlus_detcname += ".csv";
    ofstream file_muonPlus_edep_detc (edep_muonPlus_detcname, ios::out | ios::trunc);
    file_muonPlus_edep_detc<<"E(keV)"<<"\t"<<
  "eventID"<<endl;


}
////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////
    if(fSaveProcPart == true){
    G4String  Process_Particles_Name;
    Process_Particles_Name += fileName;
    Process_Particles_Name +="_Process_Particles_Name";
    Process_Particles_Name +=".csv";
    ofstream file_Process_Particles_Name (Process_Particles_Name, ios::out | ios::trunc);
    file_Process_Particles_Name<<"Edep(keV)"<<"\t"<<"Ekin(MeV)"<<"\t"<<"Proces_Name"<<"\t"<<"Particles_Name"<<"\t"<<"Parent_ID"<<"\t"<<"track_ID"<<"\t"<<
"eventID"<<endl;
    }
  /*  G4String All_electron_produc;
    All_electron_produc += fileName;
    All_electron_produc +="_All_electron_produc";
    All_electron_produc +=".csv";
    ofstream file_All_electron_produc (All_electron_produc, ios::out | ios::trunc);
    file_All_electron_produc<<"Edep(keV)"<<"\t"<<"Ekin(MeV)"<<"\t"<<"Proces_Name"<<"\t"<<"Parent_ID"<<"\t"<<"track_ID"<<"\t"<<
"eventID"<<endl;

    G4String All_gamma_produc;
    All_gamma_produc += fileName;
    All_gamma_produc +="_All_gamma_produc";
    All_gamma_produc +=".csv";
    ofstream file_All_gamma_produc (All_gamma_produc, ios::out | ios::trunc);
    file_All_gamma_produc<<"Edep(keV)"<<"\t"<<"Ekin(MeV)"<<"\t"<<"Proces_Name"<<"\t"<<"Parent_ID"<<"\t"<<"track_ID"<<"\t"<<
"eventID"<<endl;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
    G4String  Add_edep_all_electron_produc;
    Add_edep_all_electron_produc += fileName;
    Add_edep_all_electron_produc +="_Add_edep_all_electron_produc";
    Add_edep_all_electron_produc +=".csv";
    ofstream file_Add_edep_all_electron_produc (Add_edep_all_electron_produc, ios::out | ios::trunc);
    file_Add_edep_all_electron_produc<<"Edep(keV)"<<"\t"<<"Eks(MeV)"<<"\t"<<"Eke(MeV)"<<"\t"<<"count"<<"\t"<<"Ekin(MeV)"<<endl;

    G4String  Add_edep_all_gamma_produc;
    Add_edep_all_gamma_produc += fileName;
    Add_edep_all_gamma_produc +="_Add_edep_all_gamma_produc";
    Add_edep_all_gamma_produc +=".csv";
    ofstream file_Add_edep_all_gamma_produc (Add_edep_all_gamma_produc, ios::out | ios::trunc);
    file_Add_edep_all_gamma_produc<<"Edep(keV)"<<"\t"<<"Eks(MeV)"<<"\t"<<"Eke(MeV)"<<"\t"<<"count"<<"\t"<<"Ekin(MeV)"<<endl;
*/
    ////////////////////////////////////////////////////////////////////////////////////////////////
////detector
    ////////////////////////////////////////////////////////////////////////////////////////////////

    G4String  edep_all_particle_arrived_detector;
    edep_all_particle_arrived_detector += fileName;
    edep_all_particle_arrived_detector +="_edep_all_particle_arrived_detector";
    edep_all_particle_arrived_detector +=".csv";
    ofstream file_edep_all_particle_arrived_detector (edep_all_particle_arrived_detector, ios::out | ios::trunc);
    file_edep_all_particle_arrived_detector<<
              "x("<<fN_pixel_x<<")\t"<<
              "y("<<fN_pixel_y<<")\t"<<
              "Edep(keV)"<<"\t"<<
              "Ekinpost(MeV)"<<"\t"<<
              "Ekinpre(MeV)"<<"\t"<<
              "Ekin(MeV)"<<"\t"<<
              "dE/dx(keV/um)"<<"\t"<<
              "eventID"<<"\t"<<
              "Particles_Name"<<"\t"<<
              "Proces_Name"<<"\t"<<
              "Parent_ID"<<"\t"<<
              "track_ID"<<"\t"<<
              "PosX(um)"<<"\t"<<
              "PosY(um)"<<"\t"<<
              "PosZ(um)"<<"\t"<<
            endl;

    ///////////all_electron_arrived_detector////////////////////////////////////////
    G4String  edep_all_electron_arrived_detector;
    edep_all_electron_arrived_detector += fileName;
    edep_all_electron_arrived_detector +="_edep_all_electron_arrived_detector";
    edep_all_electron_arrived_detector +=".csv";
    ofstream file_edep_all_electron_arrived_detector (edep_all_electron_arrived_detector, ios::out | ios::trunc);
    file_edep_all_electron_arrived_detector<<
              "x("<<fN_pixel_x<<")\t"<<
              "y("<<fN_pixel_y<<")\t"<<
              "Edep(keV)"<<"\t"<<
              "Ekinpost(MeV)"<<"\t"<<
              "Ekinpre(MeV)"<<"\t"<<
              "Ekin(MeV)"<<"\t"<<
              "dE/dx(keV/um)"<<"\t"<<
              "eventID"<<"\t"<<

              "Proces_Name"<<"\t"<<
              "Parent_ID"<<"\t"<<
              "track_ID"<<"\t"<<

              "posX(um)\t"<<
              "posY(um)\t"<<
              "posZ(um)\t"<<
              endl;
///////////////////////////////////////////////
G4String  edep_all_anti_nu_e_arrived_detector;
edep_all_anti_nu_e_arrived_detector += fileName;
edep_all_anti_nu_e_arrived_detector +="_edep_all_anti_nu_e_arrived_detector";
edep_all_anti_nu_e_arrived_detector +=".csv";
ofstream file_edep_all_anti_nu_e_arrived_detector (edep_all_anti_nu_e_arrived_detector, ios::out | ios::trunc);
file_edep_all_anti_nu_e_arrived_detector<<
          "x("<<fN_pixel_x<<")\t"<<
          "y("<<fN_pixel_y<<")\t"<<
          "Edep(keV)"<<"\t"<<
          "Ekinpost(MeV)"<<"\t"<<
          "Ekinpre(MeV)"<<"\t"<<
          "Ekin(MeV)"<<"\t"<<
          "dE/dx(keV/um)"<<"\t"<<
          "eventID"<<"\t"<<

          "Proces_Name"<<"\t"<<
          "Parent_ID"<<"\t"<<
          "track_ID"<<"\t"<<

          "posX(um)\t"<<
          "posY(um)\t"<<
          "posZ(um)\t"<<
          endl;
///////////////////////////////////////////////

G4String  edep_all_alpha_arrived_detector;
edep_all_alpha_arrived_detector += fileName;
edep_all_alpha_arrived_detector +="_edep_all_alpha_arrived_detector";
edep_all_alpha_arrived_detector +=".csv";
ofstream file_edep_all_alpha_arrived_detector (edep_all_alpha_arrived_detector, ios::out | ios::trunc);
file_edep_all_alpha_arrived_detector<<
              "x("<<fN_pixel_x<<")\t"<<
              "y("<<fN_pixel_y<<")\t"<<
              "Edep(keV)"<<"\t"<<
              "Ekinpost(MeV)"<<"\t"<<
              "Ekinpre(MeV)"<<"\t"<<
              "Ekin(MeV)"<<"\t"<<
              "dE/dx(keV/um)"<<"\t"<<
              "eventID"<<"\t"<<

              "Proces_Name"<<"\t"<<
              "Parent_ID"<<"\t"<<
              "track_ID"<<"\t"<<

              "posX(um)\t"<<
              "posY(um)\t"<<
              "posZ(um)\t"<<
              endl;

///////
G4String  edep_all_gamma_arrived_detector;
edep_all_gamma_arrived_detector += fileName;
edep_all_gamma_arrived_detector +="_edep_all_gamma_arrived_detector";
edep_all_gamma_arrived_detector +=".csv";
ofstream file_edep_all_gamma_arrived_detector (edep_all_gamma_arrived_detector, ios::out | ios::trunc);
file_edep_all_gamma_arrived_detector<<
              "x("<<fN_pixel_x<<")\t"<<
              "y("<<fN_pixel_y<<")\t"<<
              "Edep(keV)"<<"\t"<<
              "Ekinpost(MeV)"<<"\t"<<
              "Ekinpre(MeV)"<<"\t"<<
              "Ekin(MeV)"<<"\t"<<
              "dE/dx(keV/um)"<<"\t"<<
              "eventID"<<"\t"<<

              "Proces_Name"<<"\t"<<
              "Parent_ID"<<"\t"<<
              "track_ID"<<"\t"<<

              "posX(um)\t"<<
              "posY(um)\t"<<
              "posZ(um)\t"<<
              endl;

////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

G4String  edep_zero_all_particle_arrived_detector;
edep_zero_all_particle_arrived_detector += fileName;
edep_zero_all_particle_arrived_detector +="_edep_zero_all_particle_arrived_detector";
edep_zero_all_particle_arrived_detector +=".csv";
ofstream file_edep_zero_all_particle_arrived_detector (edep_zero_all_particle_arrived_detector, ios::out | ios::trunc);
file_edep_zero_all_particle_arrived_detector<<
          "x("<<fN_pixel_x<<")\t"<<
          "y("<<fN_pixel_y<<")\t"<<
          "Edep(keV)"<<"\t"<<
          "Ekinpost(MeV)"<<"\t"<<
          "Ekinpre(MeV)"<<"\t"<<
          "Ekin(MeV)"<<"\t"<<
          "dE/dx(keV/um)"<<"\t"<<
          "eventID"<<"\t"<<
          "Particles_Name"<<"\t"<<
          "Proces_Name"<<"\t"<<
          "Parent_ID"<<"\t"<<
          "track_ID"<<"\t"<<
          "PosX(um)"<<"\t"<<
          "PosY(um)"<<"\t"<<
          "PosZ(um)"<<"\t"<<
        endl;

///////////all_electron_arrived_detector////////////////////////////////////////
G4String  edep_zero_all_electron_arrived_detector;
edep_zero_all_electron_arrived_detector += fileName;
edep_zero_all_electron_arrived_detector +="_edep_zero_all_electron_arrived_detector";
edep_zero_all_electron_arrived_detector +=".csv";
ofstream file_edep_zero_all_electron_arrived_detector (edep_zero_all_electron_arrived_detector, ios::out | ios::trunc);
file_edep_zero_all_electron_arrived_detector<<
          "x("<<fN_pixel_x<<")\t"<<
          "y("<<fN_pixel_y<<")\t"<<
          "Edep(keV)"<<"\t"<<
          "Ekinpost(MeV)"<<"\t"<<
          "Ekinpre(MeV)"<<"\t"<<
          "Ekin(MeV)"<<"\t"<<
          "dE/dx(keV/um)"<<"\t"<<
          "eventID"<<"\t"<<

          "Proces_Name"<<"\t"<<
          "Parent_ID"<<"\t"<<
          "track_ID"<<"\t"<<

          "posX(um)\t"<<
          "posY(um)\t"<<
          "posZ(um)\t"<<
          endl;
///////////////////////////////////////////////
G4String  edep_zero_all_anti_nu_e_arrived_detector;
edep_zero_all_anti_nu_e_arrived_detector += fileName;
edep_zero_all_anti_nu_e_arrived_detector +="_edep_zero_all_anti_nu_e_arrived_detector";
edep_zero_all_anti_nu_e_arrived_detector +=".csv";
ofstream file_edep_zero_all_anti_nu_e_arrived_detector (edep_zero_all_anti_nu_e_arrived_detector, ios::out | ios::trunc);
file_edep_zero_all_anti_nu_e_arrived_detector<<
      "x("<<fN_pixel_x<<")\t"<<
      "y("<<fN_pixel_y<<")\t"<<
      "Edep(keV)"<<"\t"<<
      "Ekinpost(MeV)"<<"\t"<<
      "Ekinpre(MeV)"<<"\t"<<
      "Ekin(MeV)"<<"\t"<<
      "dE/dx(keV/um)"<<"\t"<<
      "eventID"<<"\t"<<

      "Proces_Name"<<"\t"<<
      "Parent_ID"<<"\t"<<
      "track_ID"<<"\t"<<

      "posX(um)\t"<<
      "posY(um)\t"<<
      "posZ(um)\t"<<
      endl;
///////////////////////////////////////////////

G4String  edep_zero_all_alpha_arrived_detector;
edep_zero_all_alpha_arrived_detector += fileName;
edep_zero_all_alpha_arrived_detector +="_edep_zero_all_alpha_arrived_detector";
edep_zero_all_alpha_arrived_detector +=".csv";
ofstream file_edep_zero_all_alpha_arrived_detector (edep_zero_all_alpha_arrived_detector, ios::out | ios::trunc);
file_edep_zero_all_alpha_arrived_detector<<
          "x("<<fN_pixel_x<<")\t"<<
          "y("<<fN_pixel_y<<")\t"<<
          "Edep(keV)"<<"\t"<<
          "Ekinpost(MeV)"<<"\t"<<
          "Ekinpre(MeV)"<<"\t"<<
          "Ekin(MeV)"<<"\t"<<
          "dE/dx(keV/um)"<<"\t"<<
          "eventID"<<"\t"<<

          "Proces_Name"<<"\t"<<
          "Parent_ID"<<"\t"<<
          "track_ID"<<"\t"<<

          "posX(um)\t"<<
          "posY(um)\t"<<
          "posZ(um)\t"<<
          endl;

///////
G4String  edep_zero_all_gamma_arrived_detector;
edep_zero_all_gamma_arrived_detector += fileName;
edep_zero_all_gamma_arrived_detector +="_edep_zero_all_gamma_arrived_detector";
edep_zero_all_gamma_arrived_detector +=".csv";
ofstream file_edep_zero_all_gamma_arrived_detector (edep_zero_all_gamma_arrived_detector, ios::out | ios::trunc);
file_edep_zero_all_gamma_arrived_detector<<
          "x("<<fN_pixel_x<<")\t"<<
          "y("<<fN_pixel_y<<")\t"<<
          "Edep(keV)"<<"\t"<<
          "Ekinpost(MeV)"<<"\t"<<
          "Ekinpre(MeV)"<<"\t"<<
          "Ekin(MeV)"<<"\t"<<
          "dE/dx(keV/um)"<<"\t"<<
          "eventID"<<"\t"<<

          "Proces_Name"<<"\t"<<
          "Parent_ID"<<"\t"<<
          "track_ID"<<"\t"<<

          "posX(um)\t"<<
          "posY(um)\t"<<
          "posZ(um)\t"<<
          endl;


////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

G4String  edep_dif0_all_particle_arrived_detector;
edep_dif0_all_particle_arrived_detector += fileName;
edep_dif0_all_particle_arrived_detector +="_edep_dif0_all_particle_arrived_detector";
edep_dif0_all_particle_arrived_detector +=".csv";
ofstream file_edep_dif0_all_particle_arrived_detector (edep_dif0_all_particle_arrived_detector, ios::out | ios::trunc);
file_edep_dif0_all_particle_arrived_detector<<
          "x("<<fN_pixel_x<<")\t"<<
          "y("<<fN_pixel_y<<")\t"<<
          "Edep(keV)"<<"\t"<<
          "Ekinpost(MeV)"<<"\t"<<
          "Ekinpre(MeV)"<<"\t"<<
          "Ekin(MeV)"<<"\t"<<
          "dE/dx(keV/um)"<<"\t"<<
          "eventID"<<"\t"<<
          "Particles_Name"<<"\t"<<
          "Proces_Name"<<"\t"<<
          "Parent_ID"<<"\t"<<
          "track_ID"<<"\t"<<
          "PosX(um)"<<"\t"<<
          "PosY(um)"<<"\t"<<
          "PosZ(um)"<<"\t"<<
        endl;

///////////all_electron_arrived_detector////////////////////////////////////////
G4String  edep_dif0_all_electron_arrived_detector;
edep_dif0_all_electron_arrived_detector += fileName;
edep_dif0_all_electron_arrived_detector +="_edep_dif0_all_electron_arrived_detector";
edep_dif0_all_electron_arrived_detector +=".csv";
ofstream file_edep_dif0_all_electron_arrived_detector (edep_dif0_all_electron_arrived_detector, ios::out | ios::trunc);
file_edep_dif0_all_electron_arrived_detector<<
          "x("<<fN_pixel_x<<")\t"<<
          "y("<<fN_pixel_y<<")\t"<<
          "Edep(keV)"<<"\t"<<
          "Ekinpost(MeV)"<<"\t"<<
          "Ekinpre(MeV)"<<"\t"<<
          "Ekin(MeV)"<<"\t"<<
          "dE/dx(keV/um)"<<"\t"<<
          "eventID"<<"\t"<<

          "Proces_Name"<<"\t"<<
          "Parent_ID"<<"\t"<<
          "track_ID"<<"\t"<<

          "posX(um)\t"<<
          "posY(um)\t"<<
          "posZ(um)\t"<<
          endl;
///////////////////////////////////////////////
G4String  edep_dif0_all_anti_nu_e_arrived_detector;
edep_dif0_all_anti_nu_e_arrived_detector += fileName;
edep_dif0_all_anti_nu_e_arrived_detector +="_edep_dif0_all_anti_nu_e_arrived_detector";
edep_dif0_all_anti_nu_e_arrived_detector +=".csv";
ofstream file_edep_dif0_all_anti_nu_e_arrived_detector (edep_dif0_all_anti_nu_e_arrived_detector, ios::out | ios::trunc);
file_edep_dif0_all_anti_nu_e_arrived_detector<<
      "x("<<fN_pixel_x<<")\t"<<
      "y("<<fN_pixel_y<<")\t"<<
      "Edep(keV)"<<"\t"<<
      "Ekinpost(MeV)"<<"\t"<<
      "Ekinpre(MeV)"<<"\t"<<
      "Ekin(MeV)"<<"\t"<<
      "dE/dx(keV/um)"<<"\t"<<
      "eventID"<<"\t"<<

      "Proces_Name"<<"\t"<<
      "Parent_ID"<<"\t"<<
      "track_ID"<<"\t"<<

      "posX(um)\t"<<
      "posY(um)\t"<<
      "posZ(um)\t"<<
      endl;
///////////////////////////////////////////////

G4String  edep_dif0_all_alpha_arrived_detector;
edep_dif0_all_alpha_arrived_detector += fileName;
edep_dif0_all_alpha_arrived_detector +="_edep_dif0_all_alpha_arrived_detector";
edep_dif0_all_alpha_arrived_detector +=".csv";
ofstream file_edep_dif0_all_alpha_arrived_detector (edep_dif0_all_alpha_arrived_detector, ios::out | ios::trunc);
file_edep_dif0_all_alpha_arrived_detector<<
          "x("<<fN_pixel_x<<")\t"<<
          "y("<<fN_pixel_y<<")\t"<<
          "Edep(keV)"<<"\t"<<
          "Ekinpost(MeV)"<<"\t"<<
          "Ekinpre(MeV)"<<"\t"<<
          "Ekin(MeV)"<<"\t"<<
          "dE/dx(keV/um)"<<"\t"<<
          "eventID"<<"\t"<<

          "Proces_Name"<<"\t"<<
          "Parent_ID"<<"\t"<<
          "track_ID"<<"\t"<<

          "posX(um)\t"<<
          "posY(um)\t"<<
          "posZ(um)\t"<<
          endl;

///////
G4String  edep_dif0_all_gamma_arrived_detector;
edep_dif0_all_gamma_arrived_detector += fileName;
edep_dif0_all_gamma_arrived_detector +="_edep_dif0_all_gamma_arrived_detector";
edep_dif0_all_gamma_arrived_detector +=".csv";
ofstream file_edep_dif0_all_gamma_arrived_detector (edep_dif0_all_gamma_arrived_detector, ios::out | ios::trunc);
file_edep_dif0_all_gamma_arrived_detector<<
          "x("<<fN_pixel_x<<")\t"<<
          "y("<<fN_pixel_y<<")\t"<<
          "Edep(keV)"<<"\t"<<
          "Ekinpost(MeV)"<<"\t"<<
          "Ekinpre(MeV)"<<"\t"<<
          "Ekin(MeV)"<<"\t"<<
          "dE/dx(keV/um)"<<"\t"<<
          "eventID"<<"\t"<<

          "Proces_Name"<<"\t"<<
          "Parent_ID"<<"\t"<<
          "track_ID"<<"\t"<<

          "posX(um)\t"<<
          "posY(um)\t"<<
          "posZ(um)\t"<<
          endl;


////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////















    G4String  edep_electron_primary_detector;
    edep_electron_primary_detector += fileName;
    edep_electron_primary_detector +="_edep_electron_primary_detector";
    edep_electron_primary_detector +=".csv";
    ofstream file_edep_electron_primary_detector (edep_electron_primary_detector, ios::out | ios::trunc);
    file_edep_electron_primary_detector<<
              "Edep(keV)"<<"\t"<<
              "Ekin(MeV)"<<"\t"<<
              //"Proces_Name"<<"\t"<<
              "Parent_ID"<<"\t"<<
              "track_ID"<<"\t"<<
"eventID"<<endl;

    G4String  edep_electron_second_detector;
    edep_electron_second_detector += fileName;
    edep_electron_second_detector +="_edep_electron_second_detector";
    edep_electron_second_detector +=".csv";
    ofstream file_edep_electron_second_detector (edep_electron_second_detector, ios::out | ios::trunc);
    file_edep_electron_second_detector<<
              "Edep(keV)"<<"\t"<<
              "Ekin(MeV)"<<"\t"<<
              //"Proces_Name"<<"\t"<<
              "Parent_ID"<<"\t"<<
              "track_ID"<<"\t"<<
"eventID"<<endl;

    G4String  edep_gamma_second_detector;
    edep_gamma_second_detector += fileName;
    edep_gamma_second_detector +="_edep_gamma_second_detector";
    edep_gamma_second_detector +=".csv";
    ofstream file_edep_gamma_second_detector (edep_gamma_second_detector, ios::out | ios::trunc);
    file_edep_gamma_second_detector<<
              "Edep(keV)"<<"\t"<<
              "Ekin(MeV)"<<"\t"<<
              //"Proces_Name"<<"\t"<<
              "Parent_ID"<<"\t"<<
              "track_ID"<<"\t"<<
"eventID"<<endl;

if(fSource =="e" || fSource =="alpha") {
}
else{
//  G4cout << "testttttttttttttttttttttttttttt: "  << G4endl;

    ////////electron_second_parent_Sr90_detector
    G4String  edep_electron_second_parent_Sr90_detector;
    edep_electron_second_parent_Sr90_detector += fileName;
    edep_electron_second_parent_Sr90_detector +="_edep_electron_second_parent_Sr90_detector";
    edep_electron_second_parent_Sr90_detector +=".csv";
    ofstream file_edep_electron_second_parent_Sr90_detector (edep_electron_second_parent_Sr90_detector, ios::out | ios::trunc);
    file_edep_electron_second_parent_Sr90_detector<<
              "Edep(keV)"<<"\t"<<
              "Ekin(MeV)"<<"\t"<<
"eventID"<<endl;
              //"\t"<<"Proces_Name"<<endl;

    G4String  edep_electron_second_transport_parent_Sr90_detector;
    edep_electron_second_transport_parent_Sr90_detector += fileName;
    edep_electron_second_transport_parent_Sr90_detector +="_edep_electron_second_transport_parent_Sr90_detector";
    edep_electron_second_transport_parent_Sr90_detector +=".csv";
    ofstream file_edep_electron_second_transport_parent_Sr90_detector (edep_electron_second_transport_parent_Sr90_detector, ios::out | ios::trunc);
    file_edep_electron_second_transport_parent_Sr90_detector<<
              "Edep(keV)"<<"\t"<<"Ekin(MeV)"<<"\t"<<
"eventID"<<endl;

    G4String  edep_electron_second_eIoni_parent_Sr90_detector;
    edep_electron_second_eIoni_parent_Sr90_detector += fileName;
    edep_electron_second_eIoni_parent_Sr90_detector +="_edep_electron_second_eIoni_parent_Sr90_detector";
    edep_electron_second_eIoni_parent_Sr90_detector +=".csv";
    ofstream file_edep_electron_second_eIoni_parent_Sr90_detector (edep_electron_second_eIoni_parent_Sr90_detector, ios::out | ios::trunc);
    file_edep_electron_second_eIoni_parent_Sr90_detector<<
              "Edep(keV)"<<"\t"<<"Ekin(MeV)"<<"\t"<<
"eventID"<<endl;

    G4String  edep_electron_second_eBrem_parent_Sr90_detector;
    edep_electron_second_eBrem_parent_Sr90_detector += fileName;
    edep_electron_second_eBrem_parent_Sr90_detector +="_edep_electron_second_eBrem_parent_Sr90_detector";
    edep_electron_second_eBrem_parent_Sr90_detector +=".csv";
    ofstream file_edep_electron_second_eBrem_parent_Sr90_detector (edep_electron_second_eBrem_parent_Sr90_detector, ios::out | ios::trunc);
    file_edep_electron_second_eBrem_parent_Sr90_detector<<
              "Edep(keV)"<<"\t"<<"Ekin(MeV)"<<"\t"<<
"eventID"<<endl;

    G4String  edep_electron_second_msc_parent_Sr90_detector;
    edep_electron_second_msc_parent_Sr90_detector += fileName;
    edep_electron_second_msc_parent_Sr90_detector +="_edep_electron_second_msc_parent_Sr90_detector";
    edep_electron_second_msc_parent_Sr90_detector +=".csv";
    ofstream file_edep_electron_second_msc_parent_Sr90_detector (edep_electron_second_msc_parent_Sr90_detector, ios::out | ios::trunc);
    file_edep_electron_second_msc_parent_Sr90_detector<<
              "Edep(keV)"<<"\t"<<"Ekin(MeV)"<<"\t"<<
"eventID"<<endl;

////////////////////////////////////////////////////////////////////////////////
    G4String  edep_electron_second_parent_Y90_detector;
    edep_electron_second_parent_Y90_detector += fileName;
    edep_electron_second_parent_Y90_detector +="_edep_electron_second_parent_Y90_detector";
    edep_electron_second_parent_Y90_detector +=".csv";
    ofstream file_edep_electron_second_parent_Y90_detector (edep_electron_second_parent_Y90_detector, ios::out | ios::trunc);
    file_edep_electron_second_parent_Y90_detector<<
              "Edep(keV)"<<"\t"<<
              "Ekin(MeV)"<<"\t"<<
              //"Proces_Name"<<"\t"<<
              "Parent_ID"<<"\t"<<
              "track_ID"<<"\t"<<
"eventID"<<endl;

    G4String  edep_electron_second_parent_elec_detector;
    edep_electron_second_parent_elec_detector += fileName;
    edep_electron_second_parent_elec_detector +="_edep_electron_second_parent_elec_detector";
    edep_electron_second_parent_elec_detector +=".csv";
    ofstream file_edep_electron_second_parent_elec_detector (edep_electron_second_parent_elec_detector, ios::out | ios::trunc);
    file_edep_electron_second_parent_elec_detector<<
              "Edep(keV)"<<"\t"<<
              "Ekin(MeV)"<<"\t"<<
              //"Proces_Name"<<"\t"<<
              "Parent_ID"<<"\t"<<
              "track_ID"<<"\t"<<
            "eventID"<<endl;

////////////////////
}
//////////////


    ////////////////////////////////////////////////////////////////////////////////////////////////
    G4String  Add_edep_all_electron_arrived_detector;
    Add_edep_all_electron_arrived_detector += fileName;
    Add_edep_all_electron_arrived_detector +="_Add_edep_all_electron_arrived_detector";
    Add_edep_all_electron_arrived_detector +=".csv";
    ofstream file_Add_edep_all_electron_arrived_detector (Add_edep_all_electron_arrived_detector, ios::out | ios::trunc);
    file_Add_edep_all_electron_arrived_detector<<
              "Edep(keV)"<<"\t"<<
              "Eks(MeV)"<<"\t"<<
              "Eke(MeV)"<<"\t"<<
              "count"<<"\t"<<
          //      "dE/dx T"<<"\t"<<
            //    "dE/dx(keV/um)"<<"\t"<<
              "Ekin(MeV)"<<"\t"<<
            "eventID"<<endl;


//
G4String  Add_edep_all_alpha_arrived_detector;
Add_edep_all_alpha_arrived_detector += fileName;
Add_edep_all_alpha_arrived_detector +="_Add_edep_all_alpha_arrived_detector";
Add_edep_all_alpha_arrived_detector +=".csv";
ofstream file_Add_edep_all_alpha_arrived_detector (Add_edep_all_alpha_arrived_detector, ios::out | ios::trunc);
file_Add_edep_all_alpha_arrived_detector<<
          "Edep(keV)"<<"\t"<<
          "Eks(MeV)"<<"\t"<<
          "Eke(MeV)"<<"\t"<<
          "count"<<"\t"<<
      //      "dE/dx T"<<"\t"<<
        //    "dE/dx(keV/um)"<<"\t"<<
          "Ekin(MeV)"<<"\t"<<
        "eventID"<<endl;

G4String  Add_edep_all_gamma_arrived_detector;
Add_edep_all_gamma_arrived_detector += fileName;
Add_edep_all_gamma_arrived_detector +="_Add_edep_all_gamma_arrived_detector";
Add_edep_all_gamma_arrived_detector +=".csv";
ofstream file_Add_edep_all_gamma_arrived_detector (Add_edep_all_gamma_arrived_detector, ios::out | ios::trunc);
file_Add_edep_all_gamma_arrived_detector<<
          "Edep(keV)"<<"\t"<<
          "Eks(MeV)"<<"\t"<<
          "Eke(MeV)"<<"\t"<<
          "count"<<"\t"<<
      //      "dE/dx T"<<"\t"<<
        //    "dE/dx(keV/um)"<<"\t"<<
          "Ekin(MeV)"<<"\t"<<
        "eventID"<<endl;


//



    G4String  Add_edep_electron_primary_detector;
    Add_edep_electron_primary_detector += fileName;
    Add_edep_electron_primary_detector +="_Add_edep_electron_primary_detector";
    Add_edep_electron_primary_detector +=".csv";
    ofstream file_Add_edep_electron_primary_detector (Add_edep_electron_primary_detector, ios::out | ios::trunc);
    file_Add_edep_electron_primary_detector<<
              "Edep(keV)"<<"\t"<<
              "Eks(MeV)"<<"\t"<<
              "Eke(MeV)"<<"\t"<<
              "count"<<"\t"<<
              "Ekin(MeV)"<<"\t"<<
            "eventID"<<endl;


    G4String  Add_edep_electron_second_detector;
    Add_edep_electron_second_detector += fileName;
    Add_edep_electron_second_detector +="_Add_edep_electron_second_detector";
    Add_edep_electron_second_detector +=".csv";
    ofstream file_Add_edep_electron_second_detector (Add_edep_electron_second_detector, ios::out | ios::trunc);
    file_Add_edep_electron_second_detector<<
              "Edep(keV)"<<"\t"<<
              "Eks(MeV)"<<"\t"<<
              "Eke(MeV)"<<"\t"<<
              "count"<<"\t"<<
              "Ekin(MeV)"<<"\t"<<
            "eventID"<<endl;


    G4String  Add_edep_gamma_second_detector;
    Add_edep_gamma_second_detector += fileName;
    Add_edep_gamma_second_detector +="_Add_edep_gamma_second_detector";
    Add_edep_gamma_second_detector +=".csv";
    ofstream file_Add_edep_gamma_second_detector (Add_edep_gamma_second_detector, ios::out | ios::trunc);
    file_Add_edep_gamma_second_detector<<
              "Edep(keV)"<<"\t"<<
              "Eks(MeV)"<<"\t"<<
              "Eke(MeV)"<<"\t"<<
              "count"<<"\t"<<
              "Ekin(MeV)"<<"\t"<<
            "eventID"<<endl;


              G4String  Add_edep_electron_second_parent_elec_detector;
              Add_edep_electron_second_parent_elec_detector += fileName;
              Add_edep_electron_second_parent_elec_detector +="_Add_edep_electron_second_parent_elec_detector";
              Add_edep_electron_second_parent_elec_detector +=".csv";
              ofstream file_Add_edep_electron_second_parent_elec_detector (Add_edep_electron_second_parent_elec_detector, ios::out | ios::trunc);
              file_Add_edep_electron_second_parent_elec_detector<<
                        "Edep(keV)"<<"\t"<<
                        "Eks(MeV)"<<"\t"<<
                        "Eke(MeV)"<<"\t"<<
                        "count"<<"\t"<<
                        "Ekin(MeV)"<<"\t"<<
                      "eventID"<<endl;



if(fSource =="e" || fSource =="alpha") {
}
else{
//G4cout << "testttttttttttttttttttttttttttt: "  << G4endl;

    G4String  Add_edep_electron_second_parent_Sr90_detector;
    Add_edep_electron_second_parent_Sr90_detector += fileName;
    Add_edep_electron_second_parent_Sr90_detector +="_Add_edep_electron_second_parent_Sr90_detector";
    Add_edep_electron_second_parent_Sr90_detector +=".csv";
    ofstream file_Add_edep_electron_second_parent_Sr90_detector (Add_edep_electron_second_parent_Sr90_detector, ios::out | ios::trunc);
    file_Add_edep_electron_second_parent_Sr90_detector<<
              "Edep(keV)"<<"\t"<<
              "Eks(MeV)"<<"\t"<<
              "Eke(MeV)"<<"\t"<<
              "count"<<"\t"<<
              "Ekin(MeV)"<<"\t"<<
            "eventID"<<endl;


    G4String  Add_edep_electron_second_transport_parent_Sr90_detector;
    Add_edep_electron_second_transport_parent_Sr90_detector += fileName;
    Add_edep_electron_second_transport_parent_Sr90_detector +="_Add_edep_electron_second_transport_parent_Sr90_detector";
    Add_edep_electron_second_transport_parent_Sr90_detector +=".csv";
    ofstream file_Add_edep_electron_second_transport_parent_Sr90_detector (Add_edep_electron_second_transport_parent_Sr90_detector, ios::out | ios::trunc);
    file_Add_edep_electron_second_transport_parent_Sr90_detector<<
              "Edep(keV)"<<"\t"<<
              "Eks(MeV)"<<"\t"<<
              "Eke(MeV)"<<"\t"<<
              "count"<<"\t"<<
              "Ekin(MeV)"<<"\t"<<
            "eventID"<<endl;


    G4String  Add_edep_electron_second_eIoni_parent_Sr90_detector;
    Add_edep_electron_second_eIoni_parent_Sr90_detector += fileName;
    Add_edep_electron_second_eIoni_parent_Sr90_detector +="_Add_edep_electron_second_eIoni_parent_Sr90_detector";
    Add_edep_electron_second_eIoni_parent_Sr90_detector +=".csv";
    ofstream file_Add_edep_electron_second_eIoni_parent_Sr90_detector (Add_edep_electron_second_eIoni_parent_Sr90_detector, ios::out | ios::trunc);
    file_Add_edep_electron_second_eIoni_parent_Sr90_detector<<
              "Edep(keV)"<<"\t"<<
              "Eks(MeV)"<<"\t"<<
              "Eke(MeV)"<<"\t"<<
              "count"<<"\t"<<
              "Ekin(MeV)"<<"\t"<<
            "eventID"<<endl;

    G4String  Add_edep_electron_second_eBrem_parent_Sr90_detector;
    Add_edep_electron_second_eBrem_parent_Sr90_detector += fileName;
    Add_edep_electron_second_eBrem_parent_Sr90_detector +="_Add_edep_electron_second_eBrem_parent_Sr90_detector";
    Add_edep_electron_second_eBrem_parent_Sr90_detector +=".csv";
    ofstream file_Add_edep_electron_second_eBrem_parent_Sr90_detector (Add_edep_electron_second_eBrem_parent_Sr90_detector, ios::out | ios::trunc);
    file_Add_edep_electron_second_eBrem_parent_Sr90_detector<<
              "Edep(keV)"<<"\t"<<
              "Eks(MeV)"<<"\t"<<
              "Eke(MeV)"<<"\t"<<
              "count"<<"\t"<<
              "Ekin(MeV)"<<"\t"<<
            "eventID"<<endl;


    G4String  Add_edep_electron_second_msc_parent_Sr90_detector;
    Add_edep_electron_second_msc_parent_Sr90_detector += fileName;
    Add_edep_electron_second_msc_parent_Sr90_detector +="_Add_edep_electron_second_msc_parent_Sr90_detector";
    Add_edep_electron_second_msc_parent_Sr90_detector +=".csv";
    ofstream file_Add_edep_electron_second_msc_parent_Sr90_detector (Add_edep_electron_second_msc_parent_Sr90_detector, ios::out | ios::trunc);
    file_Add_edep_electron_second_msc_parent_Sr90_detector<<
              "Edep(keV)"<<"\t"<<
              "Eks(MeV)"<<"\t"<<
              "Eke(MeV)"<<"\t"<<
              "count"<<"\t"<<
              "Ekin(MeV)"<<"\t"<<
            "eventID"<<endl;


    G4String  Add_edep_electron_second_parent_Y90_detector;
    Add_edep_electron_second_parent_Y90_detector += fileName;
    Add_edep_electron_second_parent_Y90_detector +="_Add_edep_electron_second_parent_Y90_detector";
    Add_edep_electron_second_parent_Y90_detector +=".csv";
    ofstream file_Add_edep_electron_second_parent_Y90_detector (Add_edep_electron_second_parent_Y90_detector, ios::out | ios::trunc);
    file_Add_edep_electron_second_parent_Y90_detector<<
              "Edep(keV)"<<"\t"<<
              "Eks(MeV)"<<"\t"<<
              "Eke(MeV)"<<"\t"<<
              "count"<<"\t"<<
              "Ekin(MeV)"<<"\t"<<
            "eventID"<<endl;



/////////////////////////////////////
}
//////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////


    // reset accumulables to their initial values
    G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
    accumulableManager->Reset();

// Choose the Ramdom engine

    //inform the runManager to save random number seed
    //G4RunManager::GetRunManager()->SetRandomNumberStore(false);
    G4RunManager::GetRunManager()->SetRandomNumberStore(true);
    //G4RunManager::GetRunManager()->SetRandomNumberStorePerEvent(true);
    //CLHEP::HepRandom::setTheSeed(seed); G4Random::setTheSeed(seed);
    G4Random::setTheEngine(new CLHEP::RanecuEngine);

    long seeds[2];
    time_t systime = time(NULL);
    seeds[0] = (long) systime;
    seeds[1] = (long) (systime*G4UniformRand());
    G4Random::setTheSeeds(seeds);
    //CLHEP::HepRandom::setTheSeeds(seeds);
    //G4cout<< "Seed: "<< G4Random::getTheSeed()<< G4endl;verbose

//////////////////////////////////////////////////////////////////////


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* run)
{
    G4int nofEvents = run->GetNumberOfEvent();
    if (nofEvents == 0) return;

    // Merge accumulables
    G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
    accumulableManager->Merge();

    // Run conditions
    //  note: There is no primary generator action object for "master"
    //        run manager for multi-threaded mode.
    const PrimaryGeneratorAction* generatorAction
        = static_cast<const PrimaryGeneratorAction*>(
              G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
    G4String partName;
    if (generatorAction)
    {
        G4ParticleDefinition* particle
            = generatorAction->GetParticleGun()->GetParticleDefinition();
        partName = particle->GetParticleName();
    }





    // Print results
    //
    if (IsMaster())
    {
        G4cout
                << G4endl
                << "--------------------End of Global Run-----------------------"
                << G4endl
                << "  The run was "<< nofEvents<< " events ";
        //G4cout<< " * Total track length of electrons in 1st absorber: ";
        //G4cout<< fTotalTrackLength.GetValue() / mm<< " mm"<< G4endl;

    }
    else
    {
        G4cout
                << G4endl
                << "--------------------End of Local Run------------------------"
                << G4endl
                << "  The run was "<< nofEvents<< " "<< partName;
    }
    G4cout
            << "; Nb of 'good' e+ annihilations: "<< fGoodEvents.GetValue() << G4endl
            << G4endl
            << "------------------------------------------------------------"<< G4endl
            << G4endl;


/////////////////////////////////////////////////////////////////////////////////////
    // print histogram statistics
    //
    //  auto analysisManager = G4AnalysisManager::Instance();
    /*
        if ( analysisManager->GetH1(1) )
        {
            G4cout<< G4endl<< " ----> print histograms statistic ";
            if(isMaster)
            {
                G4cout<< "for the entire run "<< G4endl<< G4endl;
            }
            else
            {
                G4cout<< "for the local thread "<< G4endl<< G4endl;
            }
            G4cout<< " Edep Detec : mean = "
                  << G4BestUnit(analysisManager->GetH1(0)->mean(), "Energy")
                  << " rms = "
                  << G4BestUnit(analysisManager->GetH1(0)->rms(),  "Energy")<< G4endl;

        }
    */
    // save histograms & ntuple
    //
//    analysisManager->Write();
    //  analysisManager->CloseFile();


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RunAction::SetN_pixel_x(G4int value)
{
    fN_pixel_x = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4int RunAction::GetN_pixel_x()
{
    return fN_pixel_x;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RunAction::SetN_pixel_y(G4int value)
{
    fN_pixel_y = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4int RunAction::GetN_pixel_y()
{
    return fN_pixel_y;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::Setpixel_x(G4double value)
{
    fpixel_x = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double RunAction::Getpixel_x()
{
    return fpixel_x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RunAction::Setpixel_y(G4double value)
{
    fpixel_y = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double RunAction::Getpixel_y()
{
    return fpixel_y;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RunAction::SetSource(G4String value)
{
    fSource = value;
    //G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4String RunAction::GetSource()
{
    return fSource;
}
//...oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RunAction::SetThreshold(G4double value)
{
    fThreshold = value;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double RunAction::GetThreshold()
{
    return fThreshold;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RunAction::SetFilename(G4String value)
{
    fFilename = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4String RunAction::GetFilename()
{
    return fFilename;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RunAction::SetSaveProcePartname(G4bool value)
{
    fSaveProcPart = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool RunAction::GetSaveProcePartname()
{
    return fSaveProcPart;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
