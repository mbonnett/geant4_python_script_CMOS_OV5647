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
/// \file EventAction.cc
/// \brief Implementation of the EventAction class

#include "EventAction.hh"
#include "RunAction.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4Step.hh"

#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4THitsMap.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "EnergySD.hh"

//#include "Analysis.hh"
//#include "G4AnalysisManager.hh"
#include "Randomize.hh"
#include <iomanip>

#include <iostream>
#include <fstream>
using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(RunAction* runAction)
    : G4UserEventAction(),
      fRunAction(runAction),
      fCollID_cryst(-1),
      fenergy_SD_ID(-1),

      fAdd_Edep_all_electron_produc(0.),
      fAdd_Ekin_all_electron_produc(0.),
      fCount_all_electron_produc(0.),
      fEkin_end_all_electron_produc(0.),
      fEkin_start_all_electron_produc(0.),

      fAdd_Edep_all_gamma_produc(0.),
      fAdd_Ekin_all_gamma_produc(0.),
      fCount_all_gamma_produc(0.),
      fEkin_end_all_gamma_produc(0.),
      fEkin_start_all_gamma_produc(0.),

      fAdd_Edep_electron_second_transport_parent_Sr90(0.),
      fAdd_Ekin_electron_second_transport_parent_Sr90(0.),
      fCount_electron_second_transport_parent_Sr90(0.),
      fEkin_end_electron_second_transport_parent_Sr90(0.),
      fEkin_start_electron_second_transport_parent_Sr90(0.),

      //Detector//////////////////////////////////////////////////////////////
      fAdd_Edep_all_electron_arrived_detector(0.),
      fAdd_Ekin_all_electron_arrived_detector(0.),
      fEkin_start_all_electron_arrived_detector(0.),
      fEkin_end_all_electron_arrived_detector(0.),
      fCount_all_electron_arrived_detector(0.),

      fAdd_TrackL_all_electron_arrived_detector(0.),
      fAdd_dEdx_all_electron_arrived_detector(0.),

      fAdd_Edep_all_alpha_arrived_detector(0.),
      fAdd_Ekin_all_alpha_arrived_detector(0.),
      fEkin_start_all_alpha_arrived_detector(0.),
      fEkin_end_all_alpha_arrived_detector(0.),
      fCount_all_alpha_arrived_detector(0.),

      fAdd_TrackL_all_alpha_arrived_detector(0.),
      fAdd_dEdx_all_alpha_arrived_detector(0.),

      fAdd_Edep_all_gamma_arrived_detector(0.),
      fAdd_Ekin_all_gamma_arrived_detector(0.),
      fEkin_start_all_gamma_arrived_detector(0.),
      fEkin_end_all_gamma_arrived_detector(0.),
      fCount_all_gamma_arrived_detector(0.),

      fAdd_TrackL_all_gamma_arrived_detector(0.),
      fAdd_dEdx_all_gamma_arrived_detector(0.),




      fAdd_Edep_electron_primary_detector(0.),
      fAdd_Ekin_electron_primary_detector(0.),
      fCount_electron_primary_detector(0.),
      fEkin_end_electron_primary_detector(0.),
      fEkin_start_electron_primary_detector(0.),

      fAdd_Edep_electron_second_detector(0.),
      fAdd_Ekin_electron_second_detector(0.),
      fCount_electron_second_detector(0.),
      fEkin_end_electron_second_detector(0.),
      fEkin_start_electron_second_detector(0.),

      fAdd_Edep_gamma_second_detector(0.),
      fAdd_Ekin_gamma_second_detector(0.),
      fCount_gamma_second_detector(0.),
      fEkin_end_gamma_second_detector(0.),
      fEkin_start_gamma_second_detector(0.),

      /*  fAdd_Edep_electron_second_parent_prim_detector(0.),
      fAdd_Ekin_electron_second_parent_prim_detector(0.),
      */
      fAdd_Edep_electron_second_parent_Y90_detector(0.),
      fAdd_Ekin_electron_second_parent_Y90_detector(0.),
      fCount_electron_second_parent_Y90_detector(0.),
      fEkin_end_electron_second_parent_Y90_detector(0.),
      fEkin_start_electron_second_parent_Y90_detector(0.),

      fAdd_Edep_electron_second_parent_elec_detector(0.),
      fAdd_Ekin_electron_second_parent_elec_detector(0.),
      fCount_electron_second_parent_elec_detector(0.),
      fEkin_end_electron_second_parent_elec_detector(0.),
      fEkin_start_electron_second_parent_elec_detector(0.),

      fAdd_Edep_electron_second_parent_Sr90_detector(0.),
      fAdd_Ekin_electron_second_parent_Sr90_detector(0.),
      fCount_electron_second_parent_Sr90_detector(0.),
      fEkin_end_electron_second_parent_Sr90_detector(0.),
      fEkin_start_electron_second_parent_Sr90_detector(0.),

      fAdd_Edep_electron_second_transport_parent_Sr90_detector(0.),
      fAdd_Ekin_electron_second_transport_parent_Sr90_detector(0.),
      fCount_electron_second_transport_parent_Sr90_detector(0.),
      fEkin_end_electron_second_transport_parent_Sr90_detector(0.),
      fEkin_start_electron_second_transport_parent_Sr90_detector(0.),

      fAdd_Edep_electron_second_eIoni_parent_Sr90_detector(0.),
      fAdd_Ekin_electron_second_eIoni_parent_Sr90_detector(0.),
      fCount_electron_second_eIoni_parent_Sr90_detector(0.),
      fEkin_end_electron_second_eIoni_parent_Sr90_detector(0.),
      fEkin_start_electron_second_eIoni_parent_Sr90_detector(0.),

      fAdd_Edep_electron_second_eBrem_parent_Sr90_detector(0.),
      fAdd_Ekin_electron_second_eBrem_parent_Sr90_detector(0.),
      fCount_electron_second_eBrem_parent_Sr90_detector(0.),
      fEkin_end_electron_second_eBrem_parent_Sr90_detector(0.),
      fEkin_start_electron_second_eBrem_parent_Sr90_detector(0.),

      fAdd_Edep_electron_second_msc_parent_Sr90_detector(0.),
      fAdd_Ekin_electron_second_msc_parent_Sr90_detector(0.),
      fCount_electron_second_msc_parent_Sr90_detector(0.),
      fEkin_end_electron_second_msc_parent_Sr90_detector(0.),
      fEkin_start_electron_second_msc_parent_Sr90_detector(0.)

      /*      ///substrat
            fAdd_Edep_electron_primary_substrat(0.),
            fAdd_Edep_electron_second_substrat(0.),
            fAdd_Edep_gamma_second_substrat(0.),
            fAdd_Edep_electron_second_parent_prim_substrat(0.),
            fAdd_Edep_electron_second_parent_Y90_substrat(0.),
            fAdd_Edep_electron_second_parent_elec_substrat(0.),

            fAdd_Ekin_electron_primary_substrat(0.),
            fAdd_Ekin_electron_second_substrat(0.),
            fAdd_Ekin_gamma_second_substrat(0.),
            fAdd_Ekin_electron_second_parent_prim_substrat(0.),
            fAdd_Ekin_electron_second_parent_Y90_substrat(0.),
            fAdd_Ekin_electron_second_parent_elec_substrat(0.),



            ///World
            fAdd_Edep_electron_primary_world(0.),
            fAdd_Edep_electron_second_world(0.),
            fAdd_Edep_gamma_second_world(0.),
            fAdd_Edep_electron_second_parent_prim_world(0.),
            fAdd_Edep_electron_second_parent_Y90_world(0.),
            fAdd_Edep_electron_second_parent_elec_world(0.),

            fAdd_Ekin_electron_primary_world(0.),
            fAdd_Ekin_electron_second_world(0.),
            fAdd_Ekin_gamma_second_world(0.),
            fAdd_Ekin_electron_second_parent_prim_world(0.),
            fAdd_Ekin_electron_second_parent_Y90_world(0.),
            fAdd_Ekin_electron_second_parent_elec_world(0.)
      */

{


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{

}
///////////////////////////////////////////////////////////////////////////////////
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4THitsMap<G4double>*
EventAction::GetHitsCollection(G4int hcID,
                               const G4Event* event) const
{
    auto hitsCollection
        = static_cast<G4THitsMap<G4double>*>(
              event->GetHCofThisEvent()->GetHC(hcID));

    if ( ! hitsCollection )
    {
        G4ExceptionDescription msg;
        msg << "Cannot access hitsCollection ID " << hcID;
        G4Exception("EventAction::GetHitsCollection()",
                    "MyCode0003", FatalException, msg);
    }

    return hitsCollection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double EventAction::GetSum(G4THitsMap<G4double>* hitsMap) const
{
    G4double sumValue = 0.;
    for ( auto it : *hitsMap->GetMap() )
    {
        // hitsMap->GetMap() returns the map of std::map<G4int, G4double*>
        sumValue += *(it.second);
    }
    return sumValue;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::PrintEventStatistics(
    //G4double Detec_Edep, G4double absoTrackLength) const
    G4double Detec_Edep) const
{
    // Print event statistics
    //
    G4cout
            << "   detector: total energy_mik: "
            << std::setw(7) << G4BestUnit(Detec_Edep, "Energy")
            << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* /*evt*/)
{
    fAdd_Edep_all_electron_produc = 0.;
    fAdd_Ekin_all_electron_produc = 0.;
    fCount_all_electron_produc = 0.;
    fEkin_end_all_electron_produc = 0.;
    fEkin_start_all_electron_produc = 0.;

    fAdd_Edep_all_gamma_produc = 0.;
    fAdd_Ekin_all_gamma_produc = 0.;
    fCount_all_gamma_produc = 0.;
    fEkin_end_all_gamma_produc = 0.;
    fEkin_start_all_gamma_produc = 0.;

    fAdd_Edep_electron_second_transport_parent_Sr90 = 0.;
    fAdd_Ekin_electron_second_transport_parent_Sr90 = 0.;
    fCount_electron_second_transport_parent_Sr90 = 0.;
    fEkin_end_electron_second_transport_parent_Sr90 = 0.;
    fEkin_start_electron_second_transport_parent_Sr90 = 0.;

    //Detector//////////////////////////////////////////////////////////////
    fAdd_Edep_all_electron_arrived_detector = 0.;
    fAdd_Ekin_all_electron_arrived_detector = 0.;
    fEkin_start_all_electron_arrived_detector = 0.;
    fEkin_end_all_electron_arrived_detector = 0.;
    fCount_all_electron_arrived_detector = 0.;

    fAdd_TrackL_all_electron_arrived_detector = 0;
    fAdd_dEdx_all_electron_arrived_detector = 0;


    fAdd_Edep_all_alpha_arrived_detector = 0.;
    fAdd_Ekin_all_alpha_arrived_detector = 0.;
    fEkin_start_all_alpha_arrived_detector = 0.;
    fEkin_end_all_alpha_arrived_detector = 0.;
    fCount_all_alpha_arrived_detector = 0.;

    fAdd_TrackL_all_alpha_arrived_detector = 0.;
    fAdd_dEdx_all_alpha_arrived_detector = 0.;

    fAdd_Edep_all_gamma_arrived_detector = 0.;
    fAdd_Ekin_all_gamma_arrived_detector = 0.;
    fEkin_start_all_gamma_arrived_detector = 0.;
    fEkin_end_all_gamma_arrived_detector = 0.;
    fCount_all_gamma_arrived_detector = 0.;

    fAdd_TrackL_all_gamma_arrived_detector = 0.;
    fAdd_dEdx_all_gamma_arrived_detector = 0.;




    fAdd_Edep_electron_primary_detector = 0.;
    fAdd_Ekin_electron_primary_detector = 0.;
    fCount_electron_primary_detector = 0.;
    fEkin_end_electron_primary_detector = 0.;
    fEkin_start_electron_primary_detector = 0.;

    fAdd_Edep_electron_second_detector = 0.;
    fAdd_Ekin_electron_second_detector = 0.;
    fCount_electron_second_detector = 0.;
    fEkin_end_electron_second_detector = 0.;
    fEkin_start_electron_second_detector = 0.;

    fAdd_Edep_gamma_second_detector = 0.;
    fAdd_Ekin_gamma_second_detector = 0.;
    fCount_gamma_second_detector = 0.;
    fEkin_end_gamma_second_detector = 0.;
    fEkin_start_gamma_second_detector = 0.;

    /*  fAdd_Edep_electron_second_parent_prim_detector = 0.;
    fAdd_Ekin_electron_second_parent_prim_detector = 0.;
    */
    fAdd_Edep_electron_second_parent_Y90_detector = 0.;
    fAdd_Ekin_electron_second_parent_Y90_detector = 0.;
    fCount_electron_second_parent_Y90_detector = 0.;
    fEkin_end_electron_second_parent_Y90_detector = 0.;
    fEkin_start_electron_second_parent_Y90_detector = 0.;

    fAdd_Edep_electron_second_parent_elec_detector = 0.;
    fAdd_Ekin_electron_second_parent_elec_detector = 0.;
    fCount_electron_second_parent_elec_detector = 0.;
    fEkin_end_electron_second_parent_elec_detector = 0.;
    fEkin_start_electron_second_parent_elec_detector = 0.;

    fAdd_Edep_electron_second_parent_Sr90_detector = 0.;
    fAdd_Ekin_electron_second_parent_Sr90_detector = 0.;
    fCount_electron_second_parent_Sr90_detector = 0.;
    fEkin_end_electron_second_parent_Sr90_detector = 0.;
    fEkin_start_electron_second_parent_Sr90_detector = 0.;

    fAdd_Edep_electron_second_transport_parent_Sr90_detector = 0.;
    fAdd_Ekin_electron_second_transport_parent_Sr90_detector = 0.;
    fCount_electron_second_transport_parent_Sr90_detector = 0.;
    fEkin_end_electron_second_transport_parent_Sr90_detector = 0.;
    fEkin_start_electron_second_transport_parent_Sr90_detector = 0.;

    fAdd_Edep_electron_second_eIoni_parent_Sr90_detector = 0.;
    fAdd_Ekin_electron_second_eIoni_parent_Sr90_detector = 0.;
    fCount_electron_second_eIoni_parent_Sr90_detector = 0.;
    fEkin_end_electron_second_eIoni_parent_Sr90_detector = 0.;
    fEkin_start_electron_second_eIoni_parent_Sr90_detector = 0.;

    fAdd_Edep_electron_second_eBrem_parent_Sr90_detector = 0.;
    fAdd_Ekin_electron_second_eBrem_parent_Sr90_detector = 0.;
    fCount_electron_second_eBrem_parent_Sr90_detector = 0.;
    fEkin_end_electron_second_eBrem_parent_Sr90_detector = 0.;
    fEkin_start_electron_second_eBrem_parent_Sr90_detector = 0.;

    fAdd_Edep_electron_second_msc_parent_Sr90_detector = 0.;
    fAdd_Ekin_electron_second_msc_parent_Sr90_detector = 0.;
    fCount_electron_second_msc_parent_Sr90_detector = 0.;
    fEkin_end_electron_second_msc_parent_Sr90_detector = 0.;
    fEkin_start_electron_second_msc_parent_Sr90_detector = 0.;


    /*
        ///substrat
        fAdd_Edep_electron_primary_substrat = 0.;
        fAdd_Edep_electron_second_substrat = 0.;
        fAdd_Edep_gamma_second_substrat = 0.;
        fAdd_Edep_electron_second_parent_prim_substrat = 0.;
        fAdd_Edep_electron_second_parent_Y90_substrat = 0.;
        fAdd_Edep_electron_second_parent_elec_substrat = 0.;

        fAdd_Ekin_electron_primary_substrat = 0.;
        fAdd_Ekin_electron_second_substrat = 0.;
        fAdd_Ekin_gamma_second_substrat = 0.;
        fAdd_Ekin_electron_second_parent_prim_substrat = 0.;
        fAdd_Ekin_electron_second_parent_Y90_substrat = 0.;
        fAdd_Ekin_electron_second_parent_elec_substrat = 0.;

        ///World
        fAdd_Edep_electron_primary_world = 0.;
        fAdd_Edep_electron_second_world = 0.;
        fAdd_Edep_gamma_second_world = 0.;
        fAdd_Edep_electron_second_parent_prim_world = 0.;
        fAdd_Edep_electron_second_parent_Y90_world = 0.;
        fAdd_Edep_electron_second_parent_elec_world = 0.;

        fAdd_Ekin_electron_primary_world = 0.;
        fAdd_Ekin_electron_second_world = 0.;
        fAdd_Ekin_gamma_second_world = 0.;
        fAdd_Ekin_electron_second_parent_prim_world = 0.;
        fAdd_Ekin_electron_second_parent_Y90_world = 0.;
        fAdd_Ekin_electron_second_parent_elec_world = 0.;
      */


    feventID =  (G4EventManager::GetEventManager())->GetConstCurrentEvent()->GetEventID();

    fSource = fRunAction->GetSource();
    fFilename = fRunAction->GetFilename();

    fN_pxl_x = fRunAction->GetN_pixel_x();
    fN_pxl_y = fRunAction->GetN_pixel_y();
    fSize_pxl_x = fRunAction->Getpixel_x();
    fSize_pxl_y = fRunAction->Getpixel_y();
    fSaveProcPart = fRunAction->GetSaveProcePartname();



}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* evt )
{

//  G4RunManager::GetRunManager()->GetCurrentEvent() ;
  //G4EventManager::GetEventManager()->GetCurrentEvent();

    //G4Step* step;
    //G4double ekinpre = step->GetPreStepPoint()->GetKineticEnergy();
    //G4double ekinpost = step->GetPreStepPoint()->GetKineticEnergy();
    //  G4cout <<"energy" << ekinpost << endl;
    // Open output file

    //auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
    //auto eventID =evt->GetEventID();


    G4String   fileName;

    fileName = fRunAction->GetFilename();
    fileName += "_";
    fileName += fRunAction->GetSource();

/*
    G4String pxl_Edepemit_name;
    pxl_Edepemit_name += fileName;
    pxl_Edepemit_name +="_pixel_edep_emit.csv";
    ofstream pxl_Edepemit (pxl_Edepemit_name, ios::out | ios::trunc);
    pxl_Edepemit<<"E_dep(keV)"<<"\t"<<"E_emit(keV)"<<"\t"<<"eventID"<<endl;
*/

    ////////////////////////////////////////////////////////////////////////////////////////////////
/*
    G4String  Add_edep_all_electron_produc;
    Add_edep_all_electron_produc += fileName;
    Add_edep_all_electron_produc +="_Add_edep_all_electron_produc";
    Add_edep_all_electron_produc +=".csv";
    if(fAdd_Edep_all_electron_produc > 0 || fAdd_Ekin_all_electron_produc > 0)fN_pixel_y {fN_pixel_y
        ofstream file_Add_edep_all_electron_produc (Add_edep_all_electron_produc, ios::out | ios::app);
        if (file_Add_edep_all_electron_produc.is_open()) {
            file_Add_edep_all_electron_produc <<
                    fAdd_Edep_all_electron_produc / keV<< "\t"<<
                    fEkin_start_all_electron_produc<<"\t"<<
                    fEkin_end_all_electron_produc <<"\t"<<
                    fCount_all_electron_produc<<"\t"<<
                    fAdd_Ekin_all_electron_produc <<"\t"<<
                    feventID<<endl;

        }
        file_Add_edep_all_electron_produc.close();
    }

    G4String  Add_edep_all_gamma_produc;
    Add_edep_all_gamma_produc += fileName;
    Add_edep_all_gamma_produc +="_Add_edep_all_gamma_produc";
    Add_edep_all_gamma_produc +=".csv";
    if(fAdd_Edep_all_gamma_produc > 0 || fAdd_Ekin_all_gamma_produc > 0) {
        ofstream file_Add_edep_all_gamma_produc (Add_edep_all_gamma_produc, ios::out | ios::app);
        if (file_Add_edep_all_gamma_produc.is_open()) {
            file_Add_edep_all_gamma_produc <<
                    fAdd_Edep_all_gamma_produc / keV<< "\t"<<
                    fEkin_start_all_gamma_produc<<"\t"<<
                    fEkin_end_all_gamma_produc <<"\t"<<
                    fCount_all_gamma_produc<<"\t"<<
                    fAdd_Ekin_all_gamma_produc <<"\t"<<
                    feventID<<endl;
        }
        file_Add_edep_all_gamma_produc.close();
*/





/////Detector
    ////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
//G4cout << " Total Edep of alphas arrived in the detector: ";
//G4cout << fAdd_Edep_all_alpha_arrived_detector / keV << " MeV" << G4endl;

G4String  Add_edep_all_alpha_arrived_detector;
Add_edep_all_alpha_arrived_detector += fileName;
Add_edep_all_alpha_arrived_detector +="_Add_edep_all_alpha_arrived_detector";
Add_edep_all_alpha_arrived_detector +=".csv";

if(fAdd_Edep_all_alpha_arrived_detector > 0 || fAdd_Ekin_all_alpha_arrived_detector > 0) {
    ofstream file_Add_edep_all_alpha_arrived_detector (Add_edep_all_alpha_arrived_detector, ios::out | ios::app);
    if (file_Add_edep_all_alpha_arrived_detector.is_open()) {
        file_Add_edep_all_alpha_arrived_detector <<
                fAdd_Edep_all_alpha_arrived_detector/keV<< "\t"<<
                fEkin_start_all_alpha_arrived_detector<<"\t"<<
                fEkin_end_all_alpha_arrived_detector <<"\t"<<
                fCount_all_alpha_arrived_detector<<"\t"<<
//                    (fAdd_Edep_all_alpha_arrived_detector/keV)/fAdd_TrackL_all_alpha_arrived_detector<<"\t"<<
                //fAdd_dEdx_all_alpha_arrived_detector/keV<<"\t"<<
                fAdd_Ekin_all_alpha_arrived_detector <<"\t"<<
                feventID<<endl;
        //file_Add_edep_all_alpha_arrived_detector << fAdd_Edep_all_alpha_arrived_detector.GetValue() / keV << "\t"<< procName <<"\t"<< particleName<<"\t"<< parent_id<<"\t"<< track_id << endl;
    }
    file_Add_edep_all_alpha_arrived_detector.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////
//G4cout << " Total Edep of gammas arrived in the detector: ";
//G4cout << fAdd_Edep_all_gamma_arrived_detector / keV << " MeV" << G4endl;

G4String  Add_edep_all_gamma_arrived_detector;
Add_edep_all_gamma_arrived_detector += fileName;
Add_edep_all_gamma_arrived_detector +="_Add_edep_all_gamma_arrived_detector";
Add_edep_all_gamma_arrived_detector +=".csv";

if(fAdd_Edep_all_gamma_arrived_detector > 0 || fAdd_Ekin_all_gamma_arrived_detector > 0) {
    ofstream file_Add_edep_all_gamma_arrived_detector (Add_edep_all_gamma_arrived_detector, ios::out | ios::app);
    if (file_Add_edep_all_gamma_arrived_detector.is_open()) {
        file_Add_edep_all_gamma_arrived_detector <<
                fAdd_Edep_all_gamma_arrived_detector/keV<< "\t"<<
                fEkin_start_all_gamma_arrived_detector<<"\t"<<
                fEkin_end_all_gamma_arrived_detector <<"\t"<<
                fCount_all_gamma_arrived_detector<<"\t"<<
//                    (fAdd_Edep_all_gamma_arrived_detector/keV)/fAdd_TrackL_all_gamma_arrived_detector<<"\t"<<
                //fAdd_dEdx_all_gamma_arrived_detector/keV<<"\t"<<
                fAdd_Ekin_all_gamma_arrived_detector <<"\t"<<
                feventID<<endl;
        //file_Add_edep_all_gamma_arrived_detector << fAdd_Edep_all_gamma_arrived_detector.GetValue() / keV << "\t"<< procName <<"\t"<< particleName<<"\t"<< parent_id<<"\t"<< track_id << endl;
    }
    file_Add_edep_all_gamma_arrived_detector.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////

    //G4cout << " Total Edep of electrons primary in the detector: ";
    //G4cout << fAdd_Edep_electron_primary_detector / keV << " MeV" << G4endl;

    G4String  Add_edep_electron_primary_detector;
    Add_edep_electron_primary_detector += fileName;
    Add_edep_electron_primary_detector +="_Add_edep_electron_primary_detector";
    Add_edep_electron_primary_detector +=".csv";

    if(fAdd_Edep_electron_primary_detector > 0 || fAdd_Ekin_electron_primary_detector > 0) {
        ofstream file_Add_edep_electron_primary_detector (Add_edep_electron_primary_detector, ios::out | ios::app);
        if (file_Add_edep_electron_primary_detector.is_open()) {
            file_Add_edep_electron_primary_detector <<
                    fAdd_Edep_electron_primary_detector / keV<< "\t"<<
                    fEkin_start_electron_primary_detector<<"\t"<<
                    fEkin_end_electron_primary_detector <<"\t"<<
                    fCount_electron_primary_detector<<"\t"<<
                    fAdd_Ekin_electron_primary_detector <<"\t"<<
                    feventID<<endl;
        }
        file_Add_edep_electron_primary_detector.close();
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////
    //G4cout << " Total Edep of electrons secondary in the detector: ";
    //G4cout << fAdd_Edep_electron_second_detector / keV << " MeV" << G4endl;

    G4String  Add_edep_electron_second_detector;
    Add_edep_electron_second_detector += fileName;
    Add_edep_electron_second_detector +="_Add_edep_electron_second_detector";
    Add_edep_electron_second_detector +=".csv";

    if(fAdd_Edep_electron_second_detector > 0 || fAdd_Ekin_electron_second_detector > 0) {
        ofstream file_Add_edep_electron_second_detector (Add_edep_electron_second_detector, ios::out | ios::app);
        if (file_Add_edep_electron_second_detector.is_open()) {
            file_Add_edep_electron_second_detector <<
                    fAdd_Edep_electron_second_detector / keV<< "\t"<<
                    fEkin_start_electron_second_detector<<"\t"<<
                    fEkin_end_electron_second_detector <<"\t"<<
                    fCount_electron_second_detector<<"\t"<<
                    fAdd_Ekin_electron_second_detector <<"\t"<<
                    feventID<<endl;
        }
        file_Add_edep_electron_second_detector.close();
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////
    //G4cout << " Total Edep of gammas Secondary in the detector: ";
    //G4cout << fAdd_Edep_gamma_second_detector / keV << " MeV" << G4endl;

    G4String  Add_edep_gamma_second_detector;
    Add_edep_gamma_second_detector += fileName;
    Add_edep_gamma_second_detector +="_Add_edep_gamma_second_detector";
    Add_edep_gamma_second_detector +=".csv";

    if(fAdd_Edep_gamma_second_detector > 0 ||  fAdd_Ekin_gamma_second_detector > 0) {
        ofstream file_Add_edep_gamma_second_detector (Add_edep_gamma_second_detector, ios::out | ios::app);
        if (file_Add_edep_gamma_second_detector.is_open()) {
            file_Add_edep_gamma_second_detector <<
                    fAdd_Edep_gamma_second_detector / keV<< "\t"<<
                    fEkin_start_gamma_second_detector<<"\t"<<
                    fEkin_end_gamma_second_detector <<"\t"<<
                    fCount_gamma_second_detector<<"\t"<<
                    fAdd_Ekin_gamma_second_detector<<"\t"<<
                    feventID<<endl;
        }
        file_Add_edep_gamma_second_detector.close();
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////
    //G4cout << " Total Edep of electrons arrived in the detector: ";
    //G4cout << fAdd_Edep_all_electron_arrived_detector / keV << " MeV" << G4endl;

    G4String  Add_edep_all_electron_arrived_detector;
    Add_edep_all_electron_arrived_detector += fileName;
    Add_edep_all_electron_arrived_detector +="_Add_edep_all_electron_arrived_detector";
    Add_edep_all_electron_arrived_detector +=".csv";

    if(fAdd_Edep_all_electron_arrived_detector > 0 || fAdd_Ekin_all_electron_arrived_detector > 0) {
        ofstream file_Add_edep_all_electron_arrived_detector (Add_edep_all_electron_arrived_detector, ios::out | ios::app);
        if (file_Add_edep_all_electron_arrived_detector.is_open()) {
            file_Add_edep_all_electron_arrived_detector <<
                    fAdd_Edep_all_electron_arrived_detector/keV<< "\t"<<
                    fEkin_start_all_electron_arrived_detector<<"\t"<<
                    fEkin_end_all_electron_arrived_detector <<"\t"<<
                    fCount_all_electron_arrived_detector<<"\t"<<
//                    (fAdd_Edep_all_electron_arrived_detector/keV)/fAdd_TrackL_all_electron_arrived_detector<<"\t"<<
                    //fAdd_dEdx_all_electron_arrived_detector/keV<<"\t"<<
                    fAdd_Ekin_all_electron_arrived_detector <<"\t"<<
                    feventID<<endl;
            //file_Add_edep_all_electron_arrived_detector << fAdd_Edep_all_electron_arrived_detector.GetValue() / keV << "\t"<< procName <<"\t"<< particleName<<"\t"<< parent_id<<"\t"<< track_id << endl;
        }
        file_Add_edep_all_electron_arrived_detector.close();
    }


    ////////////////////////////////////////////////////////////////////////////////////////////////
    //G4cout << " Total Edep of electrons Secondary_but_Parent_was_a_electron in the detector: ";
    //G4cout << fAdd_Edep_electron_second_parent_elec_detector / keV << " MeV" << G4endl;

    G4String  Add_edep_electron_second_parent_elec_detector;
    Add_edep_electron_second_parent_elec_detector += fileName;
    Add_edep_electron_second_parent_elec_detector +="_Add_edep_electron_second_parent_elec_detector";
    Add_edep_electron_second_parent_elec_detector +=".csv";

    if( fAdd_Edep_electron_second_parent_elec_detector > 0 ||  fAdd_Ekin_electron_second_parent_elec_detector > 0) {
        ofstream file_Add_edep_electron_second_parent_elec_detector (Add_edep_electron_second_parent_elec_detector, ios::out | ios::app);
        if (file_Add_edep_electron_second_parent_elec_detector.is_open()) {
            file_Add_edep_electron_second_parent_elec_detector <<
                    fAdd_Edep_electron_second_parent_elec_detector / keV<< "\t"<<
                    fEkin_start_electron_second_parent_elec_detector<<"\t"<<
                    fEkin_end_electron_second_parent_elec_detector <<"\t"<<
                    fCount_electron_second_parent_elec_detector<<"\t"<<
                    fAdd_Ekin_electron_second_parent_elec_detector <<"\t"<<
                    feventID<<endl;
        }
        file_Add_edep_electron_second_parent_elec_detector.close();
    }

///////end Detector //////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////
    if(fSource =="e" || fSource =="alpha") {
    }
else{

/////////////////////////////////////////////////////////////////////////////////
  G4String Add_edep_electron_second_transport_parent_Sr90;
  Add_edep_electron_second_transport_parent_Sr90 += fileName;
  Add_edep_electron_second_transport_parent_Sr90 +="_Add_edep_electron_second_transport_parent_Sr90";
  Add_edep_electron_second_transport_parent_Sr90 +=".csv";
  if(fAdd_Edep_electron_second_transport_parent_Sr90 > 0 || fAdd_Ekin_electron_second_transport_parent_Sr90 > 0) {
      ofstream file_Add_edep_electron_second_transport_parent_Sr90 (Add_edep_electron_second_transport_parent_Sr90, ios::out | ios::app);
      if (file_Add_edep_electron_second_transport_parent_Sr90.is_open()) {
          file_Add_edep_electron_second_transport_parent_Sr90 <<
                  fAdd_Edep_electron_second_transport_parent_Sr90 / keV<< "\t"<<
                  fEkin_start_electron_second_transport_parent_Sr90<<"\t"<<
                  fEkin_end_electron_second_transport_parent_Sr90 <<"\t"<<
                  fCount_electron_second_transport_parent_Sr90<<"\t"<<
                  fAdd_Ekin_electron_second_transport_parent_Sr90 <<"\t"<<
                  feventID<<endl;
      }
      file_Add_edep_electron_second_transport_parent_Sr90.close();
  }
///////////////////////////////////////////////////////////////////////////
/////Detector
    G4String  Add_edep_electron_second_parent_Sr90_detector;
    Add_edep_electron_second_parent_Sr90_detector += fileName;
    Add_edep_electron_second_parent_Sr90_detector +="_Add_edep_electron_second_parent_Sr90_detector";
    Add_edep_electron_second_parent_Sr90_detector +=".csv";

    if( fAdd_Edep_electron_second_parent_Sr90_detector > 0 ||  fAdd_Ekin_electron_second_parent_Sr90_detector > 0) {
        ofstream file_Add_edep_electron_second_parent_Sr90_detector (Add_edep_electron_second_parent_Sr90_detector, ios::out | ios::app);
        if (file_Add_edep_electron_second_parent_Sr90_detector.is_open()) {
            file_Add_edep_electron_second_parent_Sr90_detector <<
                    fAdd_Edep_electron_second_parent_Sr90_detector / keV<< "\t"<<
                    fEkin_start_electron_second_parent_Sr90_detector<<"\t"<<
                    fEkin_end_electron_second_parent_Sr90_detector <<"\t"<<
                    fCount_electron_second_parent_Sr90_detector<<"\t"<<
                    fAdd_Ekin_electron_second_parent_Sr90_detector <<"\t"<<
                    feventID<<endl;
        }
        file_Add_edep_electron_second_parent_Sr90_detector.close();
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////
    //G4cout << " Total Edep of electrons Secondary by Transportation but Parent was Sr90 in the detector: ";
    //G4cout << fAdd_Edep_electron_second_transport_parent_Sr90_detector / keV << " MeV" << G4endl;

    G4String  Add_edep_electron_second_transport_parent_Sr90_detector;
    Add_edep_electron_second_transport_parent_Sr90_detector += fileName;
    Add_edep_electron_second_transport_parent_Sr90_detector +="_Add_edep_electron_second_transport_parent_Sr90_detector";
    Add_edep_electron_second_transport_parent_Sr90_detector +=".csv";

    if( fAdd_Edep_electron_second_transport_parent_Sr90_detector > 0 ||  fAdd_Ekin_electron_second_transport_parent_Sr90_detector > 0) {
        ofstream file_Add_edep_electron_second_transport_parent_Sr90_detector (Add_edep_electron_second_transport_parent_Sr90_detector, ios::out | ios::app);
        if (file_Add_edep_electron_second_transport_parent_Sr90_detector.is_open()) {
            file_Add_edep_electron_second_transport_parent_Sr90_detector <<
                    fAdd_Edep_electron_second_transport_parent_Sr90_detector / keV<< "\t"<<
                    fEkin_start_electron_second_transport_parent_Sr90_detector<<"\t"<<
                    fEkin_end_electron_second_transport_parent_Sr90_detector <<"\t"<<
                    fCount_electron_second_transport_parent_Sr90_detector<<"\t"<<
                    fAdd_Ekin_electron_second_transport_parent_Sr90_detector <<"\t"<<
                    feventID<<endl;
        }
        file_Add_edep_electron_second_transport_parent_Sr90_detector.close();
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////
    //G4cout << " Total Edep of electrons Secondary by eIoni but Parent was Sr90 in the detectorr: ";
    //G4cout << fAdd_Edep_electron_second_eIoni_parent_Sr90_detector / keV << " MeV" << G4endl;

    G4String  Add_edep_electron_second_eIoni_parent_Sr90_detector;
    Add_edep_electron_second_eIoni_parent_Sr90_detector += fileName;
    Add_edep_electron_second_eIoni_parent_Sr90_detector +="_Add_edep_electron_second_eIoni_parent_Sr90_detector";
    Add_edep_electron_second_eIoni_parent_Sr90_detector +=".csv";

    if( fAdd_Edep_electron_second_eIoni_parent_Sr90_detector > 0 ||  fAdd_Ekin_electron_second_eIoni_parent_Sr90_detector > 0) {
        ofstream file_Add_edep_electron_second_eIoni_parent_Sr90_detector (Add_edep_electron_second_eIoni_parent_Sr90_detector, ios::out | ios::app);
        if (file_Add_edep_electron_second_eIoni_parent_Sr90_detector.is_open()) {
            file_Add_edep_electron_second_eIoni_parent_Sr90_detector <<
                    fAdd_Edep_electron_second_eIoni_parent_Sr90_detector / keV<< "\t"<<
                    fEkin_start_electron_second_eIoni_parent_Sr90_detector<<"\t"<<
                    fEkin_end_electron_second_eIoni_parent_Sr90_detector <<"\t"<<
                    fCount_electron_second_eIoni_parent_Sr90_detector<<"\t"<<
                    fAdd_Ekin_electron_second_eIoni_parent_Sr90_detector <<"\t"<<
                    feventID<<endl;
        }
        file_Add_edep_electron_second_eIoni_parent_Sr90_detector.close();
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////
    //G4cout << " Total Edep of electrons Secondary by eBrem but Parent was Sr90 in the detector: ";
    //G4cout << fAdd_Edep_electron_second_eBrem_parent_Sr90_detector / keV << " MeV" << G4endl;

    G4String  Add_edep_electron_second_eBrem_parent_Sr90_detector;
    Add_edep_electron_second_eBrem_parent_Sr90_detector += fileName;
    Add_edep_electron_second_eBrem_parent_Sr90_detector +="_Add_edep_electron_second_eBrem_parent_Sr90_detector";
    Add_edep_electron_second_eBrem_parent_Sr90_detector +=".csv";

    if( fAdd_Edep_electron_second_eBrem_parent_Sr90_detector > 0 ||  fAdd_Ekin_electron_second_eBrem_parent_Sr90_detector > 0) {
        ofstream file_Add_edep_electron_second_eBrem_parent_Sr90_detector (Add_edep_electron_second_eBrem_parent_Sr90_detector, ios::out | ios::app);
        if (file_Add_edep_electron_second_eBrem_parent_Sr90_detector.is_open()) {
            file_Add_edep_electron_second_eBrem_parent_Sr90_detector <<
                    fAdd_Edep_electron_second_eBrem_parent_Sr90_detector / keV<< "\t"<<
                    fEkin_start_electron_second_eBrem_parent_Sr90_detector<<"\t"<<
                    fEkin_end_electron_second_eBrem_parent_Sr90_detector <<"\t"<<
                    fCount_electron_second_eBrem_parent_Sr90_detector<<"\t"<<
                    fAdd_Ekin_electron_second_eBrem_parent_Sr90_detector <<"\t"<<
                    feventID<<endl;
        }
        file_Add_edep_electron_second_eBrem_parent_Sr90_detector.close();
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////
    //G4cout << " Total Edep of electrons Secondary by msc but Parent was Sr90 in the detector: ";
    //G4cout << fAdd_Edep_electron_second_msc_parent_Sr90_detector / keV << " MeV" << G4endl;

    G4String  Add_edep_electron_second_msc_parent_Sr90_detector;
    Add_edep_electron_second_msc_parent_Sr90_detector += fileName;
    Add_edep_electron_second_msc_parent_Sr90_detector +="_Add_edep_electron_second_msc_parent_Sr90_detector";
    Add_edep_electron_second_msc_parent_Sr90_detector +=".csv";

    if( fAdd_Edep_electron_second_msc_parent_Sr90_detector > 0 ||  fAdd_Ekin_electron_second_msc_parent_Sr90_detector > 0) {
        ofstream file_Add_edep_electron_second_msc_parent_Sr90_detector (Add_edep_electron_second_msc_parent_Sr90_detector, ios::out | ios::app);
        if (file_Add_edep_electron_second_msc_parent_Sr90_detector.is_open()) {
            file_Add_edep_electron_second_msc_parent_Sr90_detector <<
                    fAdd_Edep_electron_second_msc_parent_Sr90_detector / keV<< "\t"<<
                    fEkin_start_electron_second_msc_parent_Sr90_detector<<"\t"<<
                    fEkin_end_electron_second_msc_parent_Sr90_detector <<"\t"<<
                    fCount_electron_second_msc_parent_Sr90_detector<<"\t"<<
                    fAdd_Ekin_electron_second_msc_parent_Sr90_detector <<"\t"<<
                    feventID<<endl;
        }
        file_Add_edep_electron_second_msc_parent_Sr90_detector.close();
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////
    //G4cout << " Total Edep of electrons Secondary_but_Parent_was_a_Y90 in the detector: ";
    //G4cout << fAdd_Edep_electron_second_parent_Y90_detector / keV << " MeV" << G4endl;

    G4String  Add_edep_electron_second_parent_Y90_detector;
    Add_edep_electron_second_parent_Y90_detector += fileName;
    Add_edep_electron_second_parent_Y90_detector +="_Add_edep_electron_second_parent_Y90_detector";
    Add_edep_electron_second_parent_Y90_detector +=".csv";

    if(fAdd_Edep_electron_second_parent_Y90_detector > 0 || fAdd_Ekin_electron_second_parent_Y90_detector > 0) {
        ofstream file_Add_edep_electron_second_parent_Y90_detector (Add_edep_electron_second_parent_Y90_detector, ios::out | ios::app);
        if (file_Add_edep_electron_second_parent_Y90_detector.is_open()) {
            file_Add_edep_electron_second_parent_Y90_detector <<
                    fAdd_Edep_electron_second_parent_Y90_detector / keV<< "\t"<<
                    fEkin_start_electron_second_parent_Y90_detector<<"\t"<<
                    fEkin_end_electron_second_parent_Y90_detector <<"\t"<<
                    fCount_electron_second_parent_Y90_detector<<"\t"<<
                    fAdd_Ekin_electron_second_parent_Y90_detector<<"\t"<<
                    feventID<<endl;
        }
        file_Add_edep_electron_second_parent_Y90_detector.close();
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////
}

//////////////////////////////////////////////////////////////////////////////////////////////////


//////////// get analysis manager
    //  auto analysisManager = G4AnalysisManager::Instance();

    //Hits collections
    //
    G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
    pixelSDHitsCollection * Hits = 0;
    if(!HCE) return;
    //fenergy_SD_ID  = G4SDManager::GetSDMpointer()->GetCollectionID("EnergySD/pixelSDHitsCollection");

    // Get hits collections IDs
    if (fCollID_cryst < 0)
    {
        G4SDManager* SDMan      = G4SDManager::GetSDMpointer();
        fCollID_cryst           = SDMan->GetCollectionID("crystal/edep");
        fenergy_SD_ID           = SDMan->GetCollectionID("EnergySD/pixelSDHitsCollection");

        fCollID_charg           = SDMan->GetCollectionID("crystal/chargedep");
        //fCollID_minEkin         = SDMan->GetCollectionID("crystal/edep/TrackLength");
        //fCollID_minEkin         = SDMan->GetCollectionID("crystal/minEkin");

        fCollID_TrackLength     = SDMan->GetCollectionID("crystal/TrackLength");
        fCollID_PassTrackL      = SDMan->GetCollectionID("crystal/PassTrackL");
        fCollID_nStep           = SDMan->GetCollectionID("crystal/nStep");
        //fCollID_Scond_particl   = SDMan->GetCollectionID("crystal/nScond_particl");
        //fCollID_edepScond_particl   = SDMan->GetCollectionID("crystal/edepScond_particl");
        //fCollID_minEkinScond_particl = SDMan->GetCollectionID("crystal/minEkinScond_particl");


        fCollID_alpha           = SDMan->GetCollectionID("crystal/nAlpha");
        fCollID_gamma           = SDMan->GetCollectionID("crystal/nGamma");
        fCollID_electron        = SDMan->GetCollectionID("crystal/nElectron");
        fCollID_positron        = SDMan->GetCollectionID("crystal/nPositron");
        fCollID_MuonMinus       = SDMan->GetCollectionID("crystal/nMuonMinus");
        fCollID_MuonPlus        = SDMan->GetCollectionID("crystal/nMuonPlus");

        fCollID_edepalpha       = SDMan->GetCollectionID("crystal/edepAlpha");
        fCollID_edepgamma       = SDMan->GetCollectionID("crystal/edepGamma");
        fCollID_edepelectron    = SDMan->GetCollectionID("crystal/edepElectron");
        fCollID_edeppositron    = SDMan->GetCollectionID("crystal/edepPositron");
        fCollID_edepmuonminus   = SDMan->GetCollectionID("crystal/edepMuonMinus");
        fCollID_edepmuonplus    = SDMan->GetCollectionID("crystal/edepMuonPlus");

    }





    //Hits = GetHitCollection(HCE,"pixelSDHitsCollection");
    G4double datEemited = 0.;
    if (HCE) {
        pixelSDHitsCollection* hitsCollection = static_cast<pixelSDHitsCollection*>(HCE->GetHC(fenergy_SD_ID));
        //const G4int nHits = Hits->entries();

        for (G4int j = 0; j < hitsCollection->entries(); ++j)
        //for (G4int j = 0; iHit<hitsCollection; ++j)
        {
            fdepositedEnergy = (*hitsCollection)[j]->GetDepositedEnergy();
            femittedEnergy = (*hitsCollection)[j]->GetEmittedEnergy();

            // Hacer algo con la energía depositada y detectada, por ejemplo, imprimir en la consola
            G4cout << "Deposited Energy: " << fdepositedEnergy / keV << " keV, Emitted Energy: " << femittedEnergy / MeV << " MeV" << G4endl;

            if(fdepositedEnergy == 0.0){
                G4String pxl_Edep0emit_name;
                pxl_Edep0emit_name += fileName;
                pxl_Edep0emit_name +="_pixel_edep_zero_emit.csv";
                ofstream pxl_Edep0emit (pxl_Edep0emit_name, ios::out | ios::app);
                if (pxl_Edep0emit.is_open()) {
                    pxl_Edep0emit<<fdepositedEnergy / keV<<"\t"<<femittedEnergy / keV <<"\t"<<feventID<<endl;
                }
                pxl_Edep0emit.close();
            }


            if(fdepositedEnergy >= 0.0){
                G4String pxl_Edepemit_name;
                pxl_Edepemit_name += fileName;
                pxl_Edepemit_name +="_pixel_edep_emit.csv";
                ofstream pxl_Edepemit (pxl_Edepemit_name, ios::out | ios::app);
                if (pxl_Edepemit.is_open()) {
                    pxl_Edepemit<<fdepositedEnergy / keV<<"\t"<<femittedEnergy / keV <<"\t"<<feventID<<endl;
                }
                pxl_Edepemit.close();
            }

            datEemited += femittedEnergy;
         }
    }


///////////////////////////////////////////////
// auto hc = GetHC(evt, fCollID_cryst);
    //  if ( ! hc ) return;
//////////////////////////////////////////////

    //Energy in crystals : identify 'good events'
    //
//  const G4double eThreshold = 500*keV;
    G4double eThreshold = fRunAction->GetThreshold();

    G4int nbOfFired = 0;

    G4double copyNb=-1, datEdep, datEemit;
    G4int nbx, nby;

    //G4int n_x = 2592;

    G4int n_x =fRunAction->GetN_pixel_x();

    G4THitsMap<G4double>* evtMap =
        (G4THitsMap<G4double>*)(HCE->GetHC(fCollID_cryst));

    std::map<G4int,G4double*>::iterator itr;


//G4THitsMap<G4double>* stpMap =  (G4THitsMap<G4double>*)(HCE->GetHC(fCollID_nStep));
//std::map<G4int,G4double*>::iterator itr_stp;
//itr_stp=stpMap->GetMap();

    //Store the total energy in a variable
    G4double totEdep = 0.;
    G4double totEemit = 0.;
    //std::cout << "energy emit pre" << ekinpre << std::endl;
    for (itr = evtMap->GetMap()->begin(); itr != evtMap->GetMap()->end(); itr++)
    {
        copyNb  = (itr->first);
        nby = (copyNb-1)/n_x;
        nbx =(copyNb-1)-nby*n_x;
        G4double edep = *(itr->second);
        fPlx_y = nby;
        fPlx_x = nbx;
        fEdep_plx = edep;

        //  G4double ekin = *(itr1->second);

        //Sum the energy deposited in all crystals, irrespectively of threshold.
        totEdep += edep;
        totEemit += femittedEnergy;

        if (edep >= eThreshold)
        {
            nbOfFired++; //count how many events have energy deposition above threshold

            //G4cout << "\n  pixel " << eThreshold/eV <<" : (" <<nbx <<"; " <<nby <<"): " << edep/keV << " keV ";
            datEdep= edep/keV;
            datEemit = datEemited/keV;
            //G4int id = 0;
            // write  to file
            G4String  pxledep2Dname;
            pxledep2Dname += fileName;
            pxledep2Dname +="_pixel_edep2D";
            pxledep2Dname +=".csv";
            ofstream file_pxledep2D (pxledep2Dname, ios::out | ios::app);
            if (file_pxledep2D.is_open()) {
                file_pxledep2D << nbx << "\t" << nby << "\t" << datEdep <<"\t"<<// datEemit<< "\t" <<
                            // fdepositedEnergy / keV<<"\t"<<femittedEnergy / keV <<"\t" <<
                feventID<<endl;
                // "\t"<<eventID<<endl;
                //"\t" << ekin/keV <<"\t"<<
            }
            file_pxledep2D.close();

        }
    }
    //if (nbOfFired == 2) fRunAction->CountEvent();



////////////// // Get sum values from hits collections
    //
    auto Detec_Edep = GetSum(GetHitsCollection(fCollID_cryst, evt));
    //auto Detec_minEkin = GetSum(GetHitsCollection(fCollID_minEkin, evt));

    // fill histograms
    //
    //analysisManager->FillH1(0, Detec_Edep/keV);
    //// //analysisManager->SetH1Ascii(0, true);

    // fill ntuple
    //G4int id = 1;
//   // //analysisManager->FillNtupleDColumn(id, 0, copyNb);
    //analysisManager->FillNtupleDColumn(id, 0, Detec_Edep/keV);
    G4String  edep_detcname;
    edep_detcname += fileName;
    edep_detcname += "_edep_dect";
    edep_detcname += ".csv";
    ofstream file_edep_detc (edep_detcname, ios::out | ios::app);
    if (file_edep_detc.is_open()) {
        file_edep_detc << Detec_Edep/keV <<"\t"<<  //Detec_minEkin/eV <<"\t"<<
        feventID<<endl;
    }
    file_edep_detc.close();

//    // //analysisManager->FillNtupleDColumn(id, 2, datEdep);
    // 1) total energy released in the crystals (double), MeV
    //analysisManager->FillNtupleDColumn(id, 1, totEdep/keV);
    // 2) number of crystals fired (int)
    //analysisManager->FillNtupleIColumn(id, 2, nbOfFired);
    //analysisManager>AddNtupleRow(id);



///////////////////////////////////////////////

    //print per event (modulo n)
    //
    auto eventID = evt->GetEventID();
    auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
    if ( ( printModulo > 0 ) && ( eventID % printModulo == 0 ) )
    {
        //G4cout << "---> End of event_mik: " << eventID << G4endl;
        //  PrintEventStatistics(Detec_Edep);
    }
///////////////////////////////////////////////////////////
    G4int nbxd, nbyd;
    G4double copyNbd;


////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //charge in detector
        //
        G4double charge_dep = 0.;

        evtMap = (G4THitsMap<G4double>*)(HCE->GetHC(fCollID_charg));

        for (itr = evtMap->GetMap()->begin(); itr != evtMap->GetMap()->end(); itr++)
        {
            copyNbd  = (itr->first);
            nbyd = (copyNbd-1)/n_x;
            nbxd = (copyNbd-1)-nbyd*n_x;

            charge_dep = *(itr->second);
          //  if(charge_dep!=0) {
                G4String charge_dep2Dname;
                charge_dep2Dname += fileName;
                charge_dep2Dname +="_charge_dep2D";
                charge_dep2Dname +=".csv";
                ofstream file_charge_dep2Dname (charge_dep2Dname, ios::out | ios::app);
                G4cout << "\n  Cell charge " <<": (" <<nbxd <<"; " <<nbyd <<"): " << charge_dep << " ";
                if (file_charge_dep2Dname.is_open()) {
                    file_charge_dep2Dname << nbxd << "\t" <<nbyd << "\t" << charge_dep  <<"\t"<<
                    feventID<<endl;
                }
                file_charge_dep2Dname.close();
                //    // analysisManager->FillH2(2, nbxd, nbyd, charge_dep/keV );
            //}
        }


        auto charge_Detec = GetSum(GetHitsCollection(fCollID_charg, evt));
        if(charge_Detec!=0) {
            G4String charge_dep_detcname;
            charge_dep_detcname += fileName;
            charge_dep_detcname += "_charge_dep_dect";
            charge_dep_detcname += ".csv";
            ofstream file_charge_dep_detc (charge_dep_detcname, ios::out | ios::app);
            if (file_charge_dep_detc.is_open()) {
                file_charge_dep_detc << charge_Detec <<"\t"<<
                feventID<<endl;
            }
            file_charge_dep_detc.close();
        }
/////////////////////////////////////////////////////////////////////////////////////////////////////////////


    //nGamma in detector
    //
    G4double gamma_dep = 0.;

    evtMap = (G4THitsMap<G4double>*)(HCE->GetHC(fCollID_gamma));

    for (itr = evtMap->GetMap()->begin(); itr != evtMap->GetMap()->end(); itr++)
    {
        copyNbd  = (itr->first);
        nbyd = (copyNbd-1)/n_x;
        nbxd = (copyNbd-1)-nbyd*n_x;

        gamma_dep = *(itr->second);
        //if(gamma_dep != 0) {
            //G4cout << "\n  gamma_dep " <<": (" <<nbxd <<"; " <<nbyd <<"): " << gamma_dep<< "\n ";
            G4String gamma_dep2Dname;
            gamma_dep2Dname += fileName;
            gamma_dep2Dname +="_gamma_dep2D";
            gamma_dep2Dname +=".csv";
            ofstream file_gamma_dep2Dname (gamma_dep2Dname, ios::out | ios::app);
            G4cout << "\n  gamma " <<": (" <<nbxd <<"; " <<nbyd <<"): " << gamma_dep << " ";
            if (file_gamma_dep2Dname.is_open()) {
                file_gamma_dep2Dname << nbxd << "\t" <<nbyd << "\t" << gamma_dep  <<"\t"<<
                feventID<<endl;
            }
            file_gamma_dep2Dname.close();
      //  }
    }

    auto gamma_Detec  = GetSum(GetHitsCollection(fCollID_gamma, evt));
    if(gamma_Detec  != 0) {
        G4String gamma_dep_detcname;
        gamma_dep_detcname += fileName;
        gamma_dep_detcname += "_gamma_dep_dect";
        gamma_dep_detcname += ".csv";
        ofstream file_gamma_dep_detc (gamma_dep_detcname, ios::out | ios::app);
        if (file_gamma_dep_detc.is_open()) {
            file_gamma_dep_detc << gamma_Detec  <<"\t"<<
            feventID<<endl;
        }
        file_gamma_dep_detc.close();
    }
    //edepGamma in detector
    //
    G4double edepgamma_det = 0.;

    evtMap = (G4THitsMap<G4double>*)(HCE->GetHC(fCollID_edepgamma));

    for (itr = evtMap->GetMap()->begin(); itr != evtMap->GetMap()->end(); itr++)
    {
        copyNbd  = (itr->first);
        nbyd = (copyNbd-1)/n_x;
        nbxd = (copyNbd-1)-nbyd*n_x;

        edepgamma_det = *(itr->second);
        //if(edepgamma_det>0) {
//        // analysisManager->SetH2Activation(1, true);
            G4String gamma_edep2Dname;
            gamma_edep2Dname += fileName;
            gamma_edep2Dname +="_gamma_edep2D";
            gamma_edep2Dname +=".csv";
            ofstream file_gamma_edep2D (gamma_edep2Dname, ios::out | ios::app);
            //G4cout << "\n  Edep_Gamma_det " <<": (" <<nbxd <<"; " <<nbyd <<"): " << edepgamma_det/keV << " keV ";
            if (file_gamma_edep2D.is_open()) {
                file_gamma_edep2D << nbxd << "\t" <<nbyd << "\t" << edepgamma_det/keV  <<"\t"<<
                feventID<<endl;
            }
            file_gamma_edep2D.close();
            //  // analysisManager->FillH2(1, nbxd, nbyd, edepgamma_det/keV );
        //}
    }

    //  if (edepgamma_det > 0.) fRunAction->SumDose(dose);
    auto gamma_Detec_Edep = GetSum(GetHitsCollection(fCollID_edepgamma, evt));
    if(gamma_Detec_Edep>0) {
        G4String edep_gamma_detcname;
        edep_gamma_detcname += fileName;
        edep_gamma_detcname += "_gamma_edep_dect";
        edep_gamma_detcname += ".csv";
        ofstream file_gamma_edep_detc (edep_gamma_detcname, ios::out | ios::app);
        if (file_gamma_edep_detc.is_open()) {
            file_gamma_edep_detc << gamma_Detec_Edep/keV <<"\t"<<
            feventID<<endl;
        }
        file_gamma_edep_detc.close();
    }

    //nElectron in detector
    //
    G4double electron_dep = 0.;

    evtMap = (G4THitsMap<G4double>*)(HCE->GetHC(fCollID_electron));

    for (itr = evtMap->GetMap()->begin(); itr != evtMap->GetMap()->end(); itr++)
    {
        copyNbd  = (itr->first);
        nbyd = (copyNbd-1)/n_x;
        nbxd = (copyNbd-1)-nbyd*n_x;

        electron_dep = *(itr->second);
        //if(electron_dep>0) {
            //G4cout << "\n  electron_dep " <<": (" <<nbxd <<"; " <<nbyd <<"): " << electron_dep<< "\n ";
            G4String electron_dep2Dname;
            electron_dep2Dname += fileName;
            electron_dep2Dname +="_electron_dep2D";
            electron_dep2Dname +=".csv";
            ofstream file_electron_dep2Dname (electron_dep2Dname, ios::out | ios::app);
            G4cout << "\n  elecron " <<": (" <<nbxd <<"; " <<nbyd <<"): " << electron_dep << " ";
            if (file_electron_dep2Dname.is_open()) {
                file_electron_dep2Dname << nbxd << "\t" <<nbyd << "\t" << electron_dep  <<"\t"<<
                feventID<<endl;
            }
            file_electron_dep2Dname.close();
        //}
    }

    auto electron_Detec = GetSum(GetHitsCollection(fCollID_electron, evt));
    if(electron_Detec != 0) {
        G4String electron_dep_detcname;
        electron_dep_detcname += fileName;
        electron_dep_detcname += "_electron_dep_dect";
        electron_dep_detcname += ".csv";
        ofstream file_electron_dep_detc (electron_dep_detcname, ios::out | ios::app);
        if (file_electron_dep_detc.is_open()) {
            file_electron_dep_detc << electron_Detec <<"\t"<<
            feventID<<endl;
        }
        file_electron_dep_detc.close();
    }
//////////////

    //  if (electron_dep > 0.) fRunAction->SumDose(dose);



    //edepElectron in detector
    //
    G4double edepelectron_det = 0.;

    evtMap = (G4THitsMap<G4double>*)(HCE->GetHC(fCollID_edepelectron));

    for (itr = evtMap->GetMap()->begin(); itr != evtMap->GetMap()->end(); itr++)
    {
        copyNbd  = (itr->first);
        nbyd = (copyNbd-1)/n_x;
        nbxd = (copyNbd-1)-nbyd*n_x;

        edepelectron_det = *(itr->second);

    //    G4Step* step;
      //  G4int track_id = step->GetTrack()->GetTrackID();


        if(edepelectron_det>0) {
            G4String electron_edep2Dname;
            electron_edep2Dname += fileName;
            electron_edep2Dname +="_electron_edep2D";
            electron_edep2Dname +=".csv";
            ofstream file_electron_edep2D (electron_edep2Dname, ios::out | ios::app);
            //G4cout << "\n  Edep_electr_det " <<": (" <<nbxd <<"; " <<nbyd <<"): " << edepelectron_det/keV << " keV ";
            if (file_electron_edep2D.is_open()) {
                file_electron_edep2D << nbxd << "\t" <<nbyd << "\t" << edepelectron_det/keV  <<"\t"<<
                feventID<<endl;
            }
            file_electron_edep2D.close();
            //    // analysisManager->FillH2(2, nbxd, nbyd, edepelectron_det/keV );
        }
    }
    //  if (edepelectron_det > 0.) fRunAction->SumDose(dose);
    auto electron_Detec_Edep = GetSum(GetHitsCollection(fCollID_edepelectron, evt));
    if(electron_Detec_Edep>0) {
        G4String edep_electron_detcname;
        edep_electron_detcname += fileName;
        edep_electron_detcname += "_electron_edep_dect";
        edep_electron_detcname += ".csv";
        ofstream file_electron_edep_detc (edep_electron_detcname, ios::out | ios::app);
        if (file_electron_edep_detc.is_open()) {
            file_electron_edep_detc << electron_Detec_Edep/keV<<"\t"<<
            feventID<<endl;
        }
        file_electron_edep_detc.close();
    }

    //nPositron in detector
    //
    G4double positron_det = 0.;

    evtMap = (G4THitsMap<G4double>*)(HCE->GetHC(fCollID_positron));

    for (itr = evtMap->GetMap()->begin(); itr != evtMap->GetMap()->end(); itr++)
    {
        copyNbd  = (itr->first);
        nbyd = (copyNbd-1)/n_x;
        nbxd = (copyNbd-1)-nbyd*n_x;

        positron_det = *(itr->second);
        if(positron_det>0) {
            //G4cout << "\n  positron_det " <<": (" <<nbxd <<"; " <<nbyd <<"): " << positron_det<< "\n ";
        }
    }
    //  if (positron_det > 0.) fRunAction->SumDose(dose);

    //edepPositron in detector
    //
    G4double edeppositron_det = 0.;

    evtMap = (G4THitsMap<G4double>*)(HCE->GetHC(fCollID_edeppositron));

    for (itr = evtMap->GetMap()->begin(); itr != evtMap->GetMap()->end(); itr++)
    {
        copyNbd  = (itr->first);
        nbyd = (copyNbd-1)/n_x;
        nbxd = (copyNbd-1)-nbyd*n_x;

        edeppositron_det = *(itr->second);
        if(edeppositron_det>0) {
            G4String positron_edep2Dname;
            positron_edep2Dname += fileName;
            positron_edep2Dname +="_positron_edep2D";
            positron_edep2Dname +=".csv";
            ofstream file_positron_edep2D (positron_edep2Dname, ios::out | ios::app);
            //G4cout << "\n  Edep_posistron_det " <<": (" <<nbxd <<"; " <<nbyd <<"): " << edeppositron_det/keV << " keV ";
            if (file_positron_edep2D.is_open()) {
                file_positron_edep2D << nbxd << "\t" <<nbyd << "\t" << edeppositron_det/keV  <<"\t"<<
                feventID<<endl;
            }
            file_positron_edep2D.close();
            //    // analysisManager->FillH2(3, nbxd, nbyd, edeppositron_det/keV );
        }
    }
    //  if (edeppositron_det > 0.) fRunAction->SumDose(dose);
    auto positron_Detec_Edep = GetSum(GetHitsCollection(fCollID_edeppositron, evt));
    if(positron_Detec_Edep>0) {
        G4String edep_positron_detcname;
        edep_positron_detcname += fileName;
        edep_positron_detcname += "_positron_edep_dect";
        edep_positron_detcname += ".csv";
        ofstream file_positron_edep_detc (edep_positron_detcname, ios::out | ios::app);
        if (file_positron_edep_detc.is_open()) {
            file_positron_edep_detc << positron_Detec_Edep/keV <<"\t"<<
            feventID<<endl;
        }
        file_positron_edep_detc.close();
    }

    //nAlpha in detector
    //
    G4double alpha_det = 0.;

    evtMap = (G4THitsMap<G4double>*)(HCE->GetHC(fCollID_alpha));

    for (itr = evtMap->GetMap()->begin(); itr != evtMap->GetMap()->end(); itr++)
    {
        copyNbd  = (itr->first);
        nbyd = (copyNbd-1)/n_x;
        nbxd = (copyNbd-1)-nbyd*n_x;

        if(alpha_det>0) {
            alpha_det = *(itr->second);
            //G4cout << "\n  alpha_det " <<": (" <<nbxd <<"; " <<nbyd <<"): " << alpha_det<< "\n ";
        }
    }
    //  if (alpha_det > 0.) fRunAction->SumDose(dose);

    //edepAlpha in detector
    //
    G4double edepalpha_det = 0.;

    evtMap = (G4THitsMap<G4double>*)(HCE->GetHC(fCollID_edepalpha));

    for (itr = evtMap->GetMap()->begin(); itr != evtMap->GetMap()->end(); itr++)
    {
        copyNbd  = (itr->first);
        nbyd = (copyNbd-1)/n_x;
        nbxd = (copyNbd-1)-nbyd*n_x;

        edepalpha_det = *(itr->second);
        if(edepalpha_det>0) {
            G4String alpha_edep2Dname;
            alpha_edep2Dname += fileName;
            alpha_edep2Dname +="_alpha_edep2D";
            alpha_edep2Dname +=".csv";
            ofstream file_alpha_edep2D (alpha_edep2Dname, ios::out | ios::app);
            //G4cout << "\n  Edep_alpha_det " <<": (" <<nbxd <<"; " <<nbyd <<"): " << edepalpha_det/keV << " keV ";
            if (file_alpha_edep2D.is_open()) {
                file_alpha_edep2D << nbxd << "\t" <<nbyd << "\t" << edepalpha_det/keV  <<"\t"<<
                feventID<<endl;
            }
            file_alpha_edep2D.close();
            //      // analysisManager->FillH2(4, nbxd, nbyd, edepalpha_det/keV );
        }
    }
    //  if (edepalpha_det > 0.) fRunAction->SumDose(dose);
    auto alpha_Detec_Edep = GetSum(GetHitsCollection(fCollID_edepalpha, evt));
    if(alpha_Detec_Edep>0) {
        G4String edep_alpha_detcname;
        edep_alpha_detcname += fileName;
        edep_alpha_detcname += "_alpha_edep_dect";
        edep_alpha_detcname += ".csv";
        ofstream file_alpha_edep_detc (edep_alpha_detcname, ios::out | ios::app);
        if (file_alpha_edep_detc.is_open()) {
            file_alpha_edep_detc << alpha_Detec_Edep/keV<<"\t"<<
            feventID<<endl;
        }
        file_alpha_edep_detc.close();
    }

    if(fSource=="e" || fSource=="alpha") {

    }

    else{        //nPositron in detector
        //
        G4double positron_det1 = 0.;

        evtMap = (G4THitsMap<G4double>*)(HCE->GetHC(fCollID_positron));

        for (itr = evtMap->GetMap()->begin(); itr != evtMap->GetMap()->end(); itr++)
        {
            copyNbd  = (itr->first);
            nbyd = (copyNbd-1)/n_x;
            nbxd = (copyNbd-1)-nbyd*n_x;

            positron_det1 = *(itr->second);
            if(positron_det1>0) {
                //G4cout << "\n  positron_det1 " <<": (" <<nbxd <<"; " <<nbyd <<"): " << positron_det1<< "\n ";
            }
        }
        //  if (positron_det1 > 0.) fRunAction->SumDose(dose);

        //edepPositron in detector
        //
        G4double edeppositron_det1 = 0.;

        evtMap = (G4THitsMap<G4double>*)(HCE->GetHC(fCollID_edeppositron));

        for (itr = evtMap->GetMap()->begin(); itr != evtMap->GetMap()->end(); itr++)
        {
            copyNbd  = (itr->first);
            nbyd = (copyNbd-1)/n_x;
            nbxd = (copyNbd-1)-nbyd*n_x;

            edeppositron_det1 = *(itr->second);
            if(edeppositron_det1>0) {
                G4String positron_edep2Dname;
                positron_edep2Dname += fileName;
                positron_edep2Dname +="_positron_edep2D";
                positron_edep2Dname +=".csv";
                ofstream file_positron_edep2D (positron_edep2Dname, ios::out | ios::app);
                //G4cout << "\n  Edep_posistron_det " <<": (" <<nbxd <<"; " <<nbyd <<"): " << edeppositron_det1/keV << " keV ";
                if (file_positron_edep2D.is_open()) {
                    file_positron_edep2D << nbxd << "\t" <<nbyd << "\t" << edeppositron_det1/keV  <<"\t"<<
                    feventID<<endl;
                }
                file_positron_edep2D.close();
                //    // analysisManager->FillH2(3, nbxd, nbyd, edeppositron_det1/keV );
            }
        }
        //  if (edeppositron_det1 > 0.) fRunAction->SumDose(dose);
        auto positron_Detec_Edep1 = GetSum(GetHitsCollection(fCollID_edeppositron, evt));
        if(positron_Detec_Edep1>0) {
            G4String edep_positron_detcname;
            edep_positron_detcname += fileName;
            edep_positron_detcname += "_positron_edep_dect";
            edep_positron_detcname += ".csv";
            ofstream file_positron_edep_detc (edep_positron_detcname, ios::out | ios::app);
            if (file_positron_edep_detc.is_open()) {
                file_positron_edep_detc << positron_Detec_Edep1/keV <<"\t"<<
                feventID<<endl;
            }
            file_positron_edep_detc.close();
        }
    //nMuonMinus in detector
    //
    G4double muonminus_det = 0.;

    evtMap = (G4THitsMap<G4double>*)(HCE->GetHC(fCollID_MuonMinus));

    for (itr = evtMap->GetMap()->begin(); itr != evtMap->GetMap()->end(); itr++)
    {
        copyNbd  = (itr->first);
        nbyd = (copyNbd-1)/n_x;
        nbxd = (copyNbd-1)-nbyd*n_x;

        muonminus_det = *(itr->second);
        if(muonminus_det>0) {
            //G4cout << "\n  muonMinus_det " <<": (" <<nbxd <<"; " <<nbyd <<"): " << muonminus_det<< "\n ";
        }
    }
    //  if (muonminus_det > 0.) fRunAction->SumDose(dose);

    //edepMuonminus in detector
    //
    G4double edepmuonminus_det = 0.;

    evtMap = (G4THitsMap<G4double>*)(HCE->GetHC(fCollID_edepmuonminus));

    for (itr = evtMap->GetMap()->begin(); itr != evtMap->GetMap()->end(); itr++)
    {
        copyNbd  = (itr->first);
        nbyd = (copyNbd-1)/n_x;
        nbxd = (copyNbd-1)-nbyd*n_x;

        edepmuonminus_det = *(itr->second);
        if(edepmuonminus_det>0) {
            G4String muonMinus_edep2Dname;
            muonMinus_edep2Dname += fileName;
            muonMinus_edep2Dname +="_muonMinus_edep2D";
            muonMinus_edep2Dname +=".csv";
            ofstream file_muonMinus_edep2D (muonMinus_edep2Dname, ios::out | ios::app);
            //G4cout << "\n  Edep_MuonMinus_det " <<": (" <<nbxd <<"; " <<nbyd <<"): " << edepmuonminus_det/keV << " keV ";
            if (file_muonMinus_edep2D.is_open()) {
                file_muonMinus_edep2D << nbxd << "\t" <<nbyd << "\t" << edepmuonminus_det/keV  <<"\t"<<
                feventID<<endl;
            }
            file_muonMinus_edep2D.close();
            //  // analysisManager->FillH2(5, nbxd, nbyd, edepmuonminus_det/keV );
        }
    }
    //  if (edepmuonminus_det > 0.) fRunAction->SumDose(dose);
    auto muonMinus_Detec_Edep = GetSum(GetHitsCollection(fCollID_edepmuonminus, evt));
    if(muonMinus_Detec_Edep>0) {
        G4String edep_muonMinus_detcname;
        edep_muonMinus_detcname += fileName;
        edep_muonMinus_detcname += "_muonMinus_edep_dect";
        edep_muonMinus_detcname += ".csv";
        ofstream file_muonMinus_edep_detc (edep_muonMinus_detcname, ios::out | ios::app);
        if (file_muonMinus_edep_detc.is_open()) {
            file_muonMinus_edep_detc << muonMinus_Detec_Edep/keV <<"\t"<<
            feventID<<endl;
        }
        file_muonMinus_edep_detc.close();
    }

    //nMuonPlus in detector
    //
    G4double muonplus_det = 0.;

    evtMap = (G4THitsMap<G4double>*)(HCE->GetHC(fCollID_MuonPlus));

    for (itr = evtMap->GetMap()->begin(); itr != evtMap->GetMap()->end(); itr++)
    {
        copyNbd  = (itr->first);
        nbyd = (copyNbd-1)/n_x;
        nbxd = (copyNbd-1)-nbyd*n_x;

        muonplus_det = *(itr->second);
        if(muonplus_det>0) {
            //G4cout << "\n  muonPlus_det " <<": (" <<nbxd <<"; " <<nbyd <<"): " << muonplus_det<< "\n ";
        }
    }
    //  if (muonplus_det > 0.) fRunAction->SumDose(dose);

    //edepMuonplus in detector
    //
    G4double edepmuonplus_det = 0.;

    evtMap = (G4THitsMap<G4double>*)(HCE->GetHC(fCollID_edepmuonplus));

    for (itr = evtMap->GetMap()->begin(); itr != evtMap->GetMap()->end(); itr++)
    {
        copyNbd  = (itr->first);
        nbyd = (copyNbd-1)/n_x;
        nbxd = (copyNbd-1)-nbyd*n_x;

        edepmuonplus_det = *(itr->second);
        if(edepmuonplus_det>0) {
            G4String muonPlus_edep2Dname;
            muonPlus_edep2Dname += fileName;
            muonPlus_edep2Dname +="_muonPlus_edep2D";
            muonPlus_edep2Dname +=".csv";
            ofstream file_muonPlus_edep2D (muonPlus_edep2Dname, ios::out | ios::app);
            //G4cout << "\n  Edep_MuonPlus_det " <<": (" <<nbxd <<"; " <<nbyd <<"): " << edepmuonplus_det/keV << " keV ";
            if (file_muonPlus_edep2D.is_open()) {
                file_muonPlus_edep2D << nbxd << "\t" <<nbyd << "\t" << edepmuonplus_det/keV  <<"\t"<<
                feventID<<endl;
            }
            file_muonPlus_edep2D.close();
            // // analysisManager->FillH2(6, nbxd, nbyd, edepmuonplus_det/keV );
        }
    }
    //  if (edepmuonplus_det > 0.) fRunAction->SumDose(dose);
    auto muonPlus_Detec_Edep = GetSum(GetHitsCollection(fCollID_edepmuonplus, evt));
    if(muonPlus_Detec_Edep>0) {
        G4String edep_muonPlus_detcname;
        edep_muonPlus_detcname += fileName;
        edep_muonPlus_detcname += "_muonPlus_edep_dect";
        edep_muonPlus_detcname += ".csv";
        ofstream file_muonPlus_edep_detc (edep_muonPlus_detcname, ios::out | ios::app);
        if (file_muonPlus_edep_detc.is_open()) {
            file_muonPlus_edep_detc << muonPlus_Detec_Edep/keV <<"\t"<<
            feventID<<endl;
        }
        file_muonPlus_edep_detc.close();
    }
////////////////////////////////////////////
      }
////////////////////////////////////////////////
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
