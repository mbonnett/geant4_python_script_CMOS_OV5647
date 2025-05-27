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
/// \file DetectorMessenger.cc
/// \brief Implementation of the DetectorMessenger class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
//#include "G4UIcmdWithAIntAndUnit.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4SystemOfUnits.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction* Det)
    :G4UImessenger(),
     fDetector(Det), fDetDir(0),
     fpixel_xCmd(0), fpixel_yCmd(0), fpixel_zCmd(0), fpixels_gapCmd(0),
     fN_pixel_xCmd(0), fN_pixel_yCmd(0), fDistance_SDCmd(0),fSource_typeCmd(0),
     fDetec_substraCmd(0),fSi_substrat_thickCmd(0),fDetec_shieldCmd(0),fDetc_shield_thickCmd(0),
     fDetec_SiO2_layers(0),fSiO2_layers_thinkCmd(0), fz_lensCmd(0), colr_fltrCmd(0),
     fthick_shiel_up_Cmd(0),fthick_shiel_down_Cmd(0)
{
    fDetDir = new G4UIdirectory("/pixel/det/");
    fDetDir->SetGuidance("detector construction commands");

    fpixel_xCmd =
        new G4UIcmdWithADoubleAndUnit("/pixel/det/setpixel_x", this);
    fpixel_xCmd->SetGuidance("Set pixel length X.");
    fpixel_xCmd->SetUnitCategory("Length");
    fpixel_xCmd->SetParameterName("choice",false);
    fpixel_xCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fpixel_yCmd =
        new G4UIcmdWithADoubleAndUnit("/pixel/det/setpixel_y", this);
    fpixel_yCmd->SetGuidance("Set pixel length Y.");
    fpixel_yCmd->SetUnitCategory("Length");
    fpixel_yCmd->SetParameterName("choice",false);
    fpixel_yCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fpixel_zCmd =
        new G4UIcmdWithADoubleAndUnit("/pixel/det/setpixel_z", this);
    fpixel_zCmd->SetGuidance("Set pixel length Z.");
    fpixel_zCmd->SetUnitCategory("Length");
    fpixel_zCmd->SetParameterName("choice",false);
    fpixel_zCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fpixels_gapCmd =
        new G4UIcmdWithADoubleAndUnit("/pixel/det/setpixels_gap", this);
    fpixels_gapCmd->SetGuidance("Set pixel gap.");
    fpixels_gapCmd->SetUnitCategory("Length");
    fpixels_gapCmd->SetParameterName("choice",false);
    fpixels_gapCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fN_pixel_xCmd = new G4UIcmdWithAnInteger("/pixel/det/setN_pixel_x",this);
    fN_pixel_xCmd->SetGuidance("Set number of pixels to X.");
    fN_pixel_xCmd->SetParameterName("choice",false);
    fN_pixel_xCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fN_pixel_yCmd = new G4UIcmdWithAnInteger("/pixel/det/setN_pixel_y",this);
    fN_pixel_yCmd->SetGuidance("Set number of pixels to Y.");
    fN_pixel_yCmd->SetParameterName("choice",false);
    fN_pixel_yCmd->AvailableForStates(G4State_PreInit,G4State_Idle);


    fDistance_SDCmd  =
        new G4UIcmdWithADoubleAndUnit("/pixel/det/setDistance_SD",this);
    fDistance_SDCmd ->SetGuidance("Set Distance Source to Detector.");
    fDistance_SDCmd ->SetUnitCategory("Length");
    fDistance_SDCmd->SetParameterName("choice",false);
    fDistance_SDCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

// fSource_typeCmd = new G4UIcmdWithABool("/pixel/det/setSource",this);
    fSource_typeCmd = new G4UIcmdWithAString("/pixel/det/setSource",this);
    fSource_typeCmd->SetGuidance("Select Source type.");
    fSource_typeCmd->SetParameterName("choice",false);
    fSource_typeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);


    fDetec_substraCmd = new G4UIcmdWithABool("/pixel/det/setDetec_substra",this);
    fDetec_substraCmd->SetGuidance("Active Substrate Si.");
    fDetec_substraCmd->SetParameterName("choice",false);
    fDetec_substraCmd->SetDefaultValue(false);
    fDetec_substraCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fSi_substrat_thickCmd =
        new G4UIcmdWithADoubleAndUnit("/pixel/det/setSi_substrat_thick", this);
    fSi_substrat_thickCmd->SetGuidance("Set substrat thick.");
    fSi_substrat_thickCmd->SetUnitCategory("Length");
    fSi_substrat_thickCmd->SetParameterName("choice",false);
    fSi_substrat_thickCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fDetec_shieldCmd = new G4UIcmdWithABool("/pixel/det/setDetec_shield",this);
    fDetec_shieldCmd->SetGuidance("Active Detector Shield.");
    fDetec_shieldCmd->SetParameterName("choice",false);
    fDetec_shieldCmd->SetDefaultValue(false);
    fDetec_shieldCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fDetc_shield_thickCmd =
        new G4UIcmdWithADoubleAndUnit("/pixel/det/setDetc_shield_thick", this);
    fDetc_shield_thickCmd->SetGuidance("Set shield thick.");
    fDetc_shield_thickCmd->SetUnitCategory("Length");
    fDetc_shield_thickCmd->SetParameterName("choice",false);
    fDetc_shield_thickCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

///////////////////////////////////////////////
    fDetec_SiO2_layers = new G4UIcmdWithABool("/pixel/det/setDetec_SiO2_layers",this);
    fDetec_SiO2_layers->SetGuidance("Active Detector SiO2 layers Shield.");
    fDetec_SiO2_layers->SetParameterName("choice",false);
    fDetec_SiO2_layers->SetDefaultValue(false);
    fDetec_SiO2_layers->AvailableForStates(G4State_PreInit,G4State_Idle);

    fSiO2_layers_thinkCmd =
        new G4UIcmdWithADoubleAndUnit("/pixel/det/setDetc_SiO2_layers_thick", this);
    fSiO2_layers_thinkCmd->SetGuidance("Set SiO2 layers thick.");
    fSiO2_layers_thinkCmd->SetUnitCategory("Length");
    fSiO2_layers_thinkCmd->SetParameterName("choice",false);
    fSiO2_layers_thinkCmd->AvailableForStates(G4State_PreInit,G4State_Idle);


    fz_lensCmd =
        new G4UIcmdWithADoubleAndUnit("/pixel/det/setDetc_lens_thick", this);
    fz_lensCmd->SetGuidance("Set len thick.");
    fz_lensCmd->SetUnitCategory("Length");
    fz_lensCmd->SetParameterName("choice",false);
    fz_lensCmd->AvailableForStates(G4State_PreInit,G4State_Idle);


    colr_fltrCmd =
        new G4UIcmdWithADoubleAndUnit("/pixel/det/setDetc_colrfltr_thick", this);
    colr_fltrCmd->SetGuidance("Set colr_fltr thick.");
    colr_fltrCmd->SetUnitCategory("Length");
    colr_fltrCmd->SetParameterName("choice",false);
    colr_fltrCmd->AvailableForStates(G4State_PreInit,G4State_Idle);


        fthick_shiel_down_Cmd =
            new G4UIcmdWithADoubleAndUnit("/pixel/det/set_thick_shiel_down", this);
        fthick_shiel_down_Cmd->SetGuidance("Set thick shiel down");
        fthick_shiel_down_Cmd->SetUnitCategory("Length");
        fthick_shiel_down_Cmd->SetParameterName("choice",false);
        fthick_shiel_down_Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);


        fthick_shiel_up_Cmd =
            new G4UIcmdWithADoubleAndUnit("/pixel/det/set_thick_shiel_up", this);
        fthick_shiel_up_Cmd->SetGuidance("Set thick shiel up");
        fthick_shiel_up_Cmd->SetUnitCategory("Length");
        fthick_shiel_up_Cmd->SetParameterName("choice",false);
        fthick_shiel_up_Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
    delete fpixel_xCmd;
    delete fpixel_yCmd;
    delete fpixel_zCmd;
    delete fpixels_gapCmd;

    delete fN_pixel_xCmd;
    delete fN_pixel_yCmd;

    delete fDistance_SDCmd;
    delete fSource_typeCmd;

    delete fDetec_substraCmd;
    delete fSi_substrat_thickCmd;
    delete fDetec_shieldCmd;
    delete fDetc_shield_thickCmd;

    delete fDetec_SiO2_layers;
    delete fSiO2_layers_thinkCmd;
    delete fz_lensCmd;
    delete colr_fltrCmd;

    delete fthick_shiel_down_Cmd;
    delete fthick_shiel_up_Cmd;

    delete fDetDir;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{

    if (command == fpixel_xCmd )
    {
        fDetector->Setpixel_x(fpixel_xCmd->GetNewDoubleValue(newValue));
    }

    if (command == fpixel_yCmd )
    {
        fDetector->Setpixel_y(fpixel_yCmd->GetNewDoubleValue(newValue));
    }

    if (command == fpixel_zCmd )
    {
        fDetector->Setpixel_z(fpixel_zCmd->GetNewDoubleValue(newValue));
    }

    if (command == fpixels_gapCmd )
    {
        fDetector->Setpixels_gap(fpixels_gapCmd->GetNewDoubleValue(newValue));
    }

    if (command == fN_pixel_xCmd )
    {
        fDetector->SetN_pixel_x(fN_pixel_xCmd->GetNewIntValue(newValue));
    }

    if (command == fN_pixel_yCmd )
    {
        fDetector->SetN_pixel_y(fN_pixel_yCmd->GetNewIntValue(newValue));
    }

    if (command == fDistance_SDCmd )
    {
        fDetector->SetDistance_SD(fDistance_SDCmd->GetNewDoubleValue(newValue));
    }

    if (command == fSource_typeCmd )
    {
        fDetector->SetSource(newValue);
    }

    if (command == fDetec_substraCmd )
    {
        fDetector->SetDetec_substra(fDetec_substraCmd->GetNewBoolValue(newValue));
    }

    if (command == fSi_substrat_thickCmd )
    {
        fDetector->SetSi_substrat_thick(fSi_substrat_thickCmd->GetNewDoubleValue(newValue));
    }

    if (command == fDetec_shieldCmd )
    {
        fDetector->SetDetec_shield(fDetec_shieldCmd->GetNewBoolValue(newValue));
    }

    if (command == fDetc_shield_thickCmd )
    {
        fDetector->SetDetc_shield_thick(fDetc_shield_thickCmd->GetNewDoubleValue(newValue));
    }


    if (command == fDetec_SiO2_layers )
    {
        fDetector->SetDetec_SiO2_layers(fDetec_SiO2_layers->GetNewBoolValue(newValue));
    }

    if (command == fSiO2_layers_thinkCmd )
    {
        fDetector->SetSiO2_layers_think(fSiO2_layers_thinkCmd->GetNewDoubleValue(newValue));
    }

    if (command == fz_lensCmd )
    {
        fDetector->Setz_lens(fz_lensCmd->GetNewDoubleValue(newValue));
    }

    if (command == colr_fltrCmd )
    {
        fDetector->Setz_colr_fltr(colr_fltrCmd->GetNewDoubleValue(newValue));
    }

    if (command == fthick_shiel_down_Cmd )
    {
        fDetector->Set_thick_shiel_down(fthick_shiel_down_Cmd->GetNewDoubleValue(newValue));
    }

    if (command == fthick_shiel_up_Cmd )
    {
        fDetector->Set_thick_shiel_up(fthick_shiel_up_Cmd->GetNewDoubleValue(newValue));
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
