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
/// \file EventAction.hh
/// \brief Definition of the EventAction class

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

#include "G4THitsMap.hh"


class RunAction;


/// Event action class
///
/// In EndOfEventAction() there is collected information event per event
/// from Hits Collections, and accumulated statistic for
/// RunAction::EndOfRunAction().

class EventAction : public G4UserEventAction
{
public:
    EventAction(RunAction* runAction);
    virtual ~EventAction();

    virtual void  BeginOfEventAction(const G4Event*);
    virtual void    EndOfEventAction(const G4Event*);
////////////////////////////////////////////////////////////////
    void AddE_all_electron_produc(G4double edep, G4double kenerg) {
        fAdd_Edep_all_electron_produc += edep;
        fAdd_Ekin_all_electron_produc += kenerg;

        fCount_all_electron_produc += 1;
        fEkin_end_all_electron_produc = kenerg;
        if(fCount_all_electron_produc == 1){fEkin_start_all_electron_produc = kenerg;}
    };
    void AddE_all_gamma_produc(G4double edep, G4double kenerg) {
        fAdd_Edep_all_gamma_produc += edep;
        fAdd_Ekin_all_gamma_produc += kenerg;

        fCount_all_gamma_produc += 1;
        fEkin_end_all_gamma_produc = kenerg;
        if(fCount_all_gamma_produc == 1){fEkin_start_all_gamma_produc = kenerg;}
    };
    void AddE_electron_second_transport_parent_Sr90(G4double edep, G4double kenerg) {
        fAdd_Edep_electron_second_transport_parent_Sr90 += edep;
        fAdd_Ekin_electron_second_transport_parent_Sr90 += kenerg;

        fCount_electron_second_transport_parent_Sr90 += 1;
        fEkin_end_electron_second_transport_parent_Sr90 = kenerg;
        if(fCount_electron_second_transport_parent_Sr90 == 1){fEkin_start_electron_second_transport_parent_Sr90 = kenerg;}
    };
///////////////////////////////////////////////////////////////

//Detector//////////////////////////////////////////////////////////////
    void AddE_all_electron_arrived_detector(G4double edep, G4double kenerg, G4double trackleng) {
        fAdd_Edep_all_electron_arrived_detector += edep;
        fAdd_Ekin_all_electron_arrived_detector += kenerg;
        fAdd_TrackL_all_electron_arrived_detector += trackleng;
        fAdd_dEdx_all_electron_arrived_detector += edep/trackleng;

        fCount_all_electron_arrived_detector += 1;
        fEkin_end_all_electron_arrived_detector = kenerg;
        if(fCount_all_electron_arrived_detector == 1){fEkin_start_all_electron_arrived_detector = kenerg;}
    };

    void AddE_all_anti_nu_e_arrived_detector(G4double edep, G4double kenerg, G4double trackleng) {
        fAdd_Edep_all_anti_nu_e_arrived_detector += edep;
        fAdd_Ekin_all_anti_nu_e_arrived_detector += kenerg;
        fAdd_TrackL_all_anti_nu_e_arrived_detector += trackleng;
        fAdd_dEdx_all_anti_nu_e_arrived_detector += edep/trackleng;

        fCount_all_anti_nu_e_arrived_detector += 1;
        fEkin_end_all_anti_nu_e_arrived_detector = kenerg;
        if(fCount_all_anti_nu_e_arrived_detector == 1){fEkin_start_all_anti_nu_e_arrived_detector = kenerg;}
    };


    void AddE_all_alpha_arrived_detector(G4double edep, G4double kenerg, G4double trackleng) {
        fAdd_Edep_all_alpha_arrived_detector += edep;
        fAdd_Ekin_all_alpha_arrived_detector += kenerg;
        fAdd_TrackL_all_alpha_arrived_detector += trackleng;
        fAdd_dEdx_all_alpha_arrived_detector += edep/trackleng;

        fCount_all_alpha_arrived_detector += 1;
        fEkin_end_all_alpha_arrived_detector = kenerg;
        if(fCount_all_alpha_arrived_detector == 1){fEkin_start_all_alpha_arrived_detector = kenerg;}
    };

    void AddE_all_gamma_arrived_detector(G4double edep, G4double kenerg, G4double trackleng) {
        fAdd_Edep_all_gamma_arrived_detector += edep;
        fAdd_Ekin_all_gamma_arrived_detector += kenerg;
        fAdd_TrackL_all_gamma_arrived_detector += trackleng;
        fAdd_dEdx_all_gamma_arrived_detector += edep/trackleng;

        fCount_all_gamma_arrived_detector += 1;
        fEkin_end_all_gamma_arrived_detector = kenerg;
        if(fCount_all_gamma_arrived_detector == 1){fEkin_start_all_gamma_arrived_detector = kenerg;}
    };

//Detector//////////////////////////////////////////////////////////////
    void AddE_zero_all_electron_arrived_detector(G4double edep, G4double kenerg, G4double trackleng) {
        fAdd_Edep_zero_all_electron_arrived_detector += edep;
        fAdd_Ekin_zero_all_electron_arrived_detector += kenerg;
        fAdd_TrackL_zero_all_electron_arrived_detector += trackleng;
        fAdd_dEdx_zero_all_electron_arrived_detector += edep/trackleng;

        fCount_zero_all_electron_arrived_detector += 1;
        fEkin_end_zero_all_electron_arrived_detector = kenerg;
        if(fCount_zero_all_electron_arrived_detector == 1){fEkin_start_zero_all_electron_arrived_detector = kenerg;}
    };

    void AddE_zero_all_anti_nu_e_arrived_detector(G4double edep, G4double kenerg, G4double trackleng) {
        fAdd_Edep_zero_all_anti_nu_e_arrived_detector += edep;
        fAdd_Ekin_zero_all_anti_nu_e_arrived_detector += kenerg;
        fAdd_TrackL_zero_all_anti_nu_e_arrived_detector += trackleng;
        fAdd_dEdx_zero_all_anti_nu_e_arrived_detector += edep/trackleng;

        fCount_zero_all_anti_nu_e_arrived_detector += 1;
        fEkin_end_zero_all_anti_nu_e_arrived_detector = kenerg;
        if(fCount_zero_all_anti_nu_e_arrived_detector == 1){fEkin_start_zero_all_anti_nu_e_arrived_detector = kenerg;}
    };


    void AddE_zero_all_alpha_arrived_detector(G4double edep, G4double kenerg, G4double trackleng) {
        fAdd_Edep_zero_all_alpha_arrived_detector += edep;
        fAdd_Ekin_zero_all_alpha_arrived_detector += kenerg;
        fAdd_TrackL_zero_all_alpha_arrived_detector += trackleng;
        fAdd_dEdx_zero_all_alpha_arrived_detector += edep/trackleng;

        fCount_zero_all_alpha_arrived_detector += 1;
        fEkin_end_zero_all_alpha_arrived_detector = kenerg;
        if(fCount_zero_all_alpha_arrived_detector == 1){fEkin_start_zero_all_alpha_arrived_detector = kenerg;}
    };

    void AddE_zero_all_gamma_arrived_detector(G4double edep, G4double kenerg, G4double trackleng) {
        fAdd_Edep_zero_all_gamma_arrived_detector += edep;
        fAdd_Ekin_zero_all_gamma_arrived_detector += kenerg;
        fAdd_TrackL_zero_all_gamma_arrived_detector += trackleng;
        fAdd_dEdx_zero_all_gamma_arrived_detector += edep/trackleng;

        fCount_zero_all_gamma_arrived_detector += 1;
        fEkin_end_zero_all_gamma_arrived_detector = kenerg;
        if(fCount_zero_all_gamma_arrived_detector == 1){fEkin_start_zero_all_gamma_arrived_detector = kenerg;}
    };




//Detector//////////////////////////////////////////////////////////////
    void AddE_dif0_all_electron_arrived_detector(G4double edep, G4double kenerg, G4double trackleng) {
        fAdd_Edep_dif0_all_electron_arrived_detector += edep;
        fAdd_Ekin_dif0_all_electron_arrived_detector += kenerg;
        fAdd_TrackL_dif0_all_electron_arrived_detector += trackleng;
        fAdd_dEdx_dif0_all_electron_arrived_detector += edep/trackleng;

        fCount_dif0_all_electron_arrived_detector += 1;
        fEkin_end_dif0_all_electron_arrived_detector = kenerg;
        if(fCount_dif0_all_electron_arrived_detector == 1){fEkin_start_dif0_all_electron_arrived_detector = kenerg;}
    };

    void AddE_dif0_all_anti_nu_e_arrived_detector(G4double edep, G4double kenerg, G4double trackleng) {
        fAdd_Edep_dif0_all_anti_nu_e_arrived_detector += edep;
        fAdd_Ekin_dif0_all_anti_nu_e_arrived_detector += kenerg;
        fAdd_TrackL_dif0_all_anti_nu_e_arrived_detector += trackleng;
        fAdd_dEdx_dif0_all_anti_nu_e_arrived_detector += edep/trackleng;

        fCount_dif0_all_anti_nu_e_arrived_detector += 1;
        fEkin_end_dif0_all_anti_nu_e_arrived_detector = kenerg;
        if(fCount_dif0_all_anti_nu_e_arrived_detector == 1){fEkin_start_dif0_all_anti_nu_e_arrived_detector = kenerg;}
    };

    void AddE_dif0_all_alpha_arrived_detector(G4double edep, G4double kenerg, G4double trackleng) {
        fAdd_Edep_dif0_all_alpha_arrived_detector += edep;
        fAdd_Ekin_dif0_all_alpha_arrived_detector += kenerg;
        fAdd_TrackL_dif0_all_alpha_arrived_detector += trackleng;
        fAdd_dEdx_dif0_all_alpha_arrived_detector += edep/trackleng;

        fCount_dif0_all_alpha_arrived_detector += 1;
        fEkin_end_dif0_all_alpha_arrived_detector = kenerg;
        if(fCount_dif0_all_alpha_arrived_detector == 1){fEkin_start_dif0_all_alpha_arrived_detector = kenerg;}
    };

    void AddE_dif0_all_gamma_arrived_detector(G4double edep, G4double kenerg, G4double trackleng) {
        fAdd_Edep_dif0_all_gamma_arrived_detector += edep;
        fAdd_Ekin_dif0_all_gamma_arrived_detector += kenerg;
        fAdd_TrackL_dif0_all_gamma_arrived_detector += trackleng;
        fAdd_dEdx_dif0_all_gamma_arrived_detector += edep/trackleng;

        fCount_dif0_all_gamma_arrived_detector += 1;
        fEkin_end_dif0_all_gamma_arrived_detector = kenerg;
        if(fCount_dif0_all_gamma_arrived_detector == 1){fEkin_start_dif0_all_gamma_arrived_detector = kenerg;}
    };


    void AddE_electron_primary_detector(G4double edep, G4double kenerg) {
        fAdd_Edep_electron_primary_detector += edep;
        fAdd_Ekin_electron_primary_detector += kenerg;

                fCount_electron_primary_detector += 1;
                fEkin_end_electron_primary_detector = kenerg;
                if(fCount_electron_primary_detector == 1){fEkin_start_electron_primary_detector = kenerg;}
    };
    void AddE_electron_second_detector(G4double edep, G4double kenerg) {
        fAdd_Edep_electron_second_detector += edep;
        fAdd_Ekin_electron_second_detector += kenerg;

                fCount_electron_second_detector += 1;
                fEkin_end_electron_second_detector = kenerg;
                if(fCount_electron_second_detector == 1){fEkin_start_electron_second_detector = kenerg;}
    };
    void AddE_gamma_second_detector(G4double edep, G4double kenerg) {
        fAdd_Edep_gamma_second_detector += edep;
        fAdd_Ekin_gamma_second_detector += kenerg;

                fCount_gamma_second_detector += 1;
                fEkin_end_gamma_second_detector = kenerg;
                if(fCount_gamma_second_detector == 1){fEkin_start_gamma_second_detector = kenerg;}
    };
    void AddE_electron_second_parent_Y90_detector(G4double edep, G4double kenerg) {
        fAdd_Edep_electron_second_parent_Y90_detector += edep;
        fAdd_Ekin_electron_second_parent_Y90_detector += kenerg;

                fCount_electron_second_parent_Y90_detector += 1;
                fEkin_end_electron_second_parent_Y90_detector = kenerg;
                if(fCount_electron_second_parent_Y90_detector == 1){fEkin_start_electron_second_parent_Y90_detector = kenerg;}
    };
    void AddE_electron_second_parent_elec_detector(G4double edep, G4double kenerg) {
        fAdd_Edep_electron_second_parent_elec_detector += edep;
        fAdd_Ekin_electron_second_parent_elec_detector += kenerg;

                fCount_electron_second_parent_elec_detector += 1;
                fEkin_end_electron_second_parent_elec_detector = kenerg;
                if(fCount_electron_second_parent_elec_detector == 1){fEkin_start_electron_second_parent_elec_detector = kenerg;}
    };
    void AddE_electron_second_parent_Sr90_detector(G4double edep, G4double kenerg) {
        fAdd_Edep_electron_second_parent_Sr90_detector += edep;
        fAdd_Ekin_electron_second_parent_Sr90_detector += kenerg;

                fCount_electron_second_parent_Sr90_detector += 1;
                fEkin_end_electron_second_parent_Sr90_detector = kenerg;
                if(fCount_electron_second_parent_Sr90_detector == 1){fEkin_start_electron_second_parent_Sr90_detector = kenerg;}
    };
    void AddE_electron_second_transport_parent_Sr90_detector(G4double edep, G4double kenerg) {
        fAdd_Edep_electron_second_transport_parent_Sr90_detector += edep;
        fAdd_Ekin_electron_second_transport_parent_Sr90_detector += kenerg;

                fCount_electron_second_transport_parent_Sr90_detector += 1;
                fEkin_end_electron_second_transport_parent_Sr90_detector = kenerg;
                if(fCount_electron_second_transport_parent_Sr90_detector == 1){fEkin_start_electron_second_transport_parent_Sr90_detector = kenerg;}
    };
    void AddE_electron_second_eIoni_parent_Sr90_detector(G4double edep, G4double kenerg) {
        fAdd_Edep_electron_second_eIoni_parent_Sr90_detector += edep;
        fAdd_Ekin_electron_second_eIoni_parent_Sr90_detector += kenerg;

                fCount_electron_second_eIoni_parent_Sr90_detector += 1;
                fEkin_end_electron_second_eIoni_parent_Sr90_detector = kenerg;
                if(fCount_electron_second_eIoni_parent_Sr90_detector == 1){fEkin_start_electron_second_eIoni_parent_Sr90_detector = kenerg;}
    };
    void AddE_electron_second_eBrem_parent_Sr90_detector(G4double edep, G4double kenerg) {
        fAdd_Edep_electron_second_eBrem_parent_Sr90_detector += edep;
        fAdd_Ekin_electron_second_eBrem_parent_Sr90_detector += kenerg;

                fCount_electron_second_eBrem_parent_Sr90_detector += 1;
                fEkin_end_electron_second_eBrem_parent_Sr90_detector = kenerg;
                if(fCount_electron_second_eBrem_parent_Sr90_detector == 1){fEkin_start_electron_second_eBrem_parent_Sr90_detector = kenerg;}
    };
    void AddE_electron_second_msc_parent_Sr90_detector(G4double edep, G4double kenerg) {
        fAdd_Edep_electron_second_msc_parent_Sr90_detector += edep;
        fAdd_Ekin_electron_second_msc_parent_Sr90_detector += kenerg;

                fCount_electron_second_msc_parent_Sr90_detector += 1;
                fEkin_end_electron_second_msc_parent_Sr90_detector = kenerg;
                if(fCount_electron_second_msc_parent_Sr90_detector == 1){fEkin_start_electron_second_msc_parent_Sr90_detector = kenerg;}
    };


public:
    G4int feventID;

    G4String    fSource;
    G4String    fFilename;
    G4int fPlx_y;
    G4int fPlx_x;
    G4double fEdep_plx;

    G4int fN_pxl_x;
    G4int fN_pxl_y;
    G4double fSize_pxl_x;
    G4double fSize_pxl_y;

    G4bool	fSaveProcPart;

    G4double fdepositedEnergy;
    G4double femittedEnergy;


private:

    // methods
    G4THitsMap<G4double>* GetHitsCollection(G4int hcID,
                                            const G4Event* event) const;
    G4double GetSum(G4THitsMap<G4double>* hitsMap) const;
    void PrintEventStatistics(G4double Detec_Edep) const;

    G4double GetTotal(const G4THitsMap<G4double> &map) const;
    G4double FindMinimum(const G4THitsMap<G4double> &map) const;

private:
    // data members

	RunAction*  fRunAction;
    G4int fCollID_cryst;
    G4int fenergy_SD_ID;
    G4int fCollID_charg;

    G4int fCollID_minEkin;

    G4int fCollID_TrackLength;
    G4int fCollID_PassTrackL;
    G4int fCollID_nStep;

    G4int fCollID_Scond_particl;
    G4int fCollID_edepScond_particl;
    G4int fCollID_minEkinScond_particl;

    G4int fCollID_alpha;
    G4int fCollID_gamma;
    G4int fCollID_electron;
    G4int fCollID_positron;
    G4int fCollID_MuonMinus;
    G4int fCollID_MuonPlus;

    G4int fCollID_edepalpha;
    G4int fCollID_edepgamma;
    G4int fCollID_edepelectron;
    G4int fCollID_edeppositron;
    G4int fCollID_edepmuonminus;
    G4int fCollID_edepmuonplus;

    G4int fCollID_minEkinalpha;
    G4int fCollID_minEkingamma;
    G4int fCollID_minEkinelectron;
    G4int fCollID_minEkinpositron;
    G4int fCollID_minEkinmuonminus;
    G4int fCollID_minEkinmuonplus;


private:

    G4double fAdd_Edep_all_electron_produc;
    G4double fAdd_Ekin_all_electron_produc;
    G4double fCount_all_electron_produc;
    G4double fEkin_end_all_electron_produc;
    G4double fEkin_start_all_electron_produc;

    G4double fAdd_Edep_all_gamma_produc;
    G4double fAdd_Ekin_all_gamma_produc;
    G4double fCount_all_gamma_produc;
    G4double fEkin_end_all_gamma_produc;
    G4double fEkin_start_all_gamma_produc;

    G4double fAdd_Edep_electron_second_transport_parent_Sr90;
    G4double fAdd_Ekin_electron_second_transport_parent_Sr90;
    G4double fCount_electron_second_transport_parent_Sr90;
    G4double fEkin_end_electron_second_transport_parent_Sr90;
    G4double fEkin_start_electron_second_transport_parent_Sr90;

//Detector//////////////////////////////////////////////////////////////
    G4double fAdd_Edep_all_electron_arrived_detector;
    G4double fAdd_Ekin_all_electron_arrived_detector;
    G4double fEkin_start_all_electron_arrived_detector;
    G4double fEkin_end_all_electron_arrived_detector;
    G4double fCount_all_electron_arrived_detector;

    G4double fAdd_TrackL_all_electron_arrived_detector;
    G4double fAdd_dEdx_all_electron_arrived_detector;

    G4double fAdd_Edep_all_anti_nu_e_arrived_detector;
    G4double fAdd_Ekin_all_anti_nu_e_arrived_detector;
    G4double fEkin_start_all_anti_nu_e_arrived_detector;
    G4double fEkin_end_all_anti_nu_e_arrived_detector;
    G4double fCount_all_anti_nu_e_arrived_detector;

    G4double fAdd_TrackL_all_anti_nu_e_arrived_detector;
    G4double fAdd_dEdx_all_anti_nu_e_arrived_detector;



    G4double fAdd_Edep_all_alpha_arrived_detector;
    G4double fAdd_Ekin_all_alpha_arrived_detector;
    G4double fEkin_start_all_alpha_arrived_detector;
    G4double fEkin_end_all_alpha_arrived_detector;
    G4double fCount_all_alpha_arrived_detector;

    G4double fAdd_TrackL_all_alpha_arrived_detector;
    G4double fAdd_dEdx_all_alpha_arrived_detector;

    G4double fAdd_Edep_all_gamma_arrived_detector;
    G4double fAdd_Ekin_all_gamma_arrived_detector;
    G4double fEkin_start_all_gamma_arrived_detector;
    G4double fEkin_end_all_gamma_arrived_detector;
    G4double fCount_all_gamma_arrived_detector;

    G4double fAdd_TrackL_all_gamma_arrived_detector;
    G4double fAdd_dEdx_all_gamma_arrived_detector;

///
//Detector//////////////////////////////////////////////////////////////
    G4double fAdd_Edep_zero_all_electron_arrived_detector;
    G4double fAdd_Ekin_zero_all_electron_arrived_detector;
    G4double fEkin_start_zero_all_electron_arrived_detector;
    G4double fEkin_end_zero_all_electron_arrived_detector;
    G4double fCount_zero_all_electron_arrived_detector;

    G4double fAdd_TrackL_zero_all_electron_arrived_detector;
    G4double fAdd_dEdx_zero_all_electron_arrived_detector;

    G4double fAdd_Edep_zero_all_anti_nu_e_arrived_detector;
    G4double fAdd_Ekin_zero_all_anti_nu_e_arrived_detector;
    G4double fEkin_start_zero_all_anti_nu_e_arrived_detector;
    G4double fEkin_end_zero_all_anti_nu_e_arrived_detector;
    G4double fCount_zero_all_anti_nu_e_arrived_detector;

    G4double fAdd_TrackL_zero_all_anti_nu_e_arrived_detector;
    G4double fAdd_dEdx_zero_all_anti_nu_e_arrived_detector;



    G4double fAdd_Edep_zero_all_alpha_arrived_detector;
    G4double fAdd_Ekin_zero_all_alpha_arrived_detector;
    G4double fEkin_start_zero_all_alpha_arrived_detector;
    G4double fEkin_end_zero_all_alpha_arrived_detector;
    G4double fCount_zero_all_alpha_arrived_detector;

    G4double fAdd_TrackL_zero_all_alpha_arrived_detector;
    G4double fAdd_dEdx_zero_all_alpha_arrived_detector;

    G4double fAdd_Edep_zero_all_gamma_arrived_detector;
    G4double fAdd_Ekin_zero_all_gamma_arrived_detector;
    G4double fEkin_start_zero_all_gamma_arrived_detector;
    G4double fEkin_end_zero_all_gamma_arrived_detector;
    G4double fCount_zero_all_gamma_arrived_detector;

    G4double fAdd_TrackL_zero_all_gamma_arrived_detector;
    G4double fAdd_dEdx_zero_all_gamma_arrived_detector;

////

//Detector//////////////////////////////////////////////////////////////
    G4double fAdd_Edep_dif0_all_electron_arrived_detector;
    G4double fAdd_Ekin_dif0_all_electron_arrived_detector;
    G4double fEkin_start_dif0_all_electron_arrived_detector;
    G4double fEkin_end_dif0_all_electron_arrived_detector;
    G4double fCount_dif0_all_electron_arrived_detector;

    G4double fAdd_TrackL_dif0_all_electron_arrived_detector;
    G4double fAdd_dEdx_dif0_all_electron_arrived_detector;

    G4double fAdd_Edep_dif0_all_anti_nu_e_arrived_detector;
    G4double fAdd_Ekin_dif0_all_anti_nu_e_arrived_detector;
    G4double fEkin_start_dif0_all_anti_nu_e_arrived_detector;
    G4double fEkin_end_dif0_all_anti_nu_e_arrived_detector;
    G4double fCount_dif0_all_anti_nu_e_arrived_detector;

    G4double fAdd_TrackL_dif0_all_anti_nu_e_arrived_detector;
    G4double fAdd_dEdx_dif0_all_anti_nu_e_arrived_detector;



    G4double fAdd_Edep_dif0_all_alpha_arrived_detector;
    G4double fAdd_Ekin_dif0_all_alpha_arrived_detector;
    G4double fEkin_start_dif0_all_alpha_arrived_detector;
    G4double fEkin_end_dif0_all_alpha_arrived_detector;
    G4double fCount_dif0_all_alpha_arrived_detector;

    G4double fAdd_TrackL_dif0_all_alpha_arrived_detector;
    G4double fAdd_dEdx_dif0_all_alpha_arrived_detector;

    G4double fAdd_Edep_dif0_all_gamma_arrived_detector;
    G4double fAdd_Ekin_dif0_all_gamma_arrived_detector;
    G4double fEkin_start_dif0_all_gamma_arrived_detector;
    G4double fEkin_end_dif0_all_gamma_arrived_detector;
    G4double fCount_dif0_all_gamma_arrived_detector;

    G4double fAdd_TrackL_dif0_all_gamma_arrived_detector;
    G4double fAdd_dEdx_dif0_all_gamma_arrived_detector;



    G4double fAdd_Edep_electron_primary_detector;
    G4double fAdd_Ekin_electron_primary_detector;
    G4double fCount_electron_primary_detector;
    G4double fEkin_end_electron_primary_detector;
    G4double fEkin_start_electron_primary_detector;

    G4double fAdd_Edep_electron_second_detector;
    G4double fAdd_Ekin_electron_second_detector;
    G4double fCount_electron_second_detector;
    G4double fEkin_end_electron_second_detector;
    G4double fEkin_start_electron_second_detector;

    G4double fAdd_Edep_gamma_second_detector;
    G4double fAdd_Ekin_gamma_second_detector;
    G4double fCount_gamma_second_detector;
    G4double fEkin_end_gamma_second_detector;
    G4double fEkin_start_gamma_second_detector;

    G4double fAdd_Edep_electron_second_parent_Y90_detector;
    G4double fAdd_Ekin_electron_second_parent_Y90_detector;
    G4double fCount_electron_second_parent_Y90_detector;
    G4double fEkin_end_electron_second_parent_Y90_detector;
    G4double fEkin_start_electron_second_parent_Y90_detector;

    G4double fAdd_Edep_electron_second_parent_elec_detector;
    G4double fAdd_Ekin_electron_second_parent_elec_detector;
    G4double fCount_electron_second_parent_elec_detector;
    G4double fEkin_end_electron_second_parent_elec_detector;
    G4double fEkin_start_electron_second_parent_elec_detector;

    G4double fAdd_Edep_electron_second_parent_Sr90_detector;
    G4double fAdd_Ekin_electron_second_parent_Sr90_detector;
    G4double fCount_electron_second_parent_Sr90_detector;
    G4double fEkin_end_electron_second_parent_Sr90_detector;
    G4double fEkin_start_electron_second_parent_Sr90_detector;

    G4double fAdd_Edep_electron_second_transport_parent_Sr90_detector;
    G4double fAdd_Ekin_electron_second_transport_parent_Sr90_detector;
    G4double fCount_electron_second_transport_parent_Sr90_detector;
    G4double fEkin_end_electron_second_transport_parent_Sr90_detector;
    G4double fEkin_start_electron_second_transport_parent_Sr90_detector;

    G4double fAdd_Edep_electron_second_eIoni_parent_Sr90_detector;
    G4double fAdd_Ekin_electron_second_eIoni_parent_Sr90_detector;
    G4double fCount_electron_second_eIoni_parent_Sr90_detector;
    G4double fEkin_end_electron_second_eIoni_parent_Sr90_detector;
    G4double fEkin_start_electron_second_eIoni_parent_Sr90_detector;

    G4double fAdd_Edep_electron_second_eBrem_parent_Sr90_detector;
    G4double fAdd_Ekin_electron_second_eBrem_parent_Sr90_detector;
    G4double fCount_electron_second_eBrem_parent_Sr90_detector;
    G4double fEkin_end_electron_second_eBrem_parent_Sr90_detector;
    G4double fEkin_start_electron_second_eBrem_parent_Sr90_detector;

    G4double fAdd_Edep_electron_second_msc_parent_Sr90_detector;
    G4double fAdd_Ekin_electron_second_msc_parent_Sr90_detector;
    G4double fCount_electron_second_msc_parent_Sr90_detector;
    G4double fEkin_end_electron_second_msc_parent_Sr90_detector;
    G4double fEkin_start_electron_second_msc_parent_Sr90_detector;

 ///////////////////////////////////////////////////////////////////////
     G4double fEdep_electron_arrived_detector;
     G4double fEkin_electron_arrived_detector;


/*
//Substrat//////////////////////////////////////////////////////////////
    G4double fAdd_Edep_electron_primary_substrat;
    G4double fAdd_Ekin_electron_primary_substrat;

    G4double fAdd_Edep_electron_second_substrat;
    G4double fAdd_Ekin_electron_second_substrat;

    G4double fAdd_Edep_gamma_second_substrat;
    G4double fAdd_Ekin_gamma_second_substrat;

    G4double fAdd_Edep_electron_second_parent_prim_substrat;
    G4double fAdd_Ekin_electron_second_parent_prim_substrat;

    G4double fAdd_Edep_electron_second_parent_Y90_substrat;
    G4double fAdd_Ekin_electron_second_parent_Y90_substrat;

    G4double fAdd_Edep_electron_second_parent_elec_substrat;
    G4double fAdd_Ekin_electron_second_parent_elec_substrat;

//World//////////////////////////////////////////////////////////////
    G4double fAdd_Edep_electron_primary_world;
    G4double fAdd_Ekin_electron_primary_world;

    G4double fAdd_Edep_electron_second_world;
    G4double fAdd_Ekin_electron_second_world;

    G4double fAdd_Edep_gamma_second_world;
    G4double fAdd_Ekin_gamma_second_world;

    G4double fAdd_Edep_electron_second_parent_prim_world;
    G4double fAdd_Ekin_electron_second_parent_prim_world;

    G4double fAdd_Edep_electron_second_parent_Y90_world;
    G4double fAdd_Ekin_electron_second_parent_Y90_world;

    G4double fAdd_Edep_electron_second_parent_elec_world;
    G4double fAdd_Ekin_electron_second_parent_elec_world;
*/
////////////////////////////////////////////////////////////////////////


};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
