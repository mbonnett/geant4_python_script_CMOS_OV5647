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
/// \file RunActionMessenger.cc
/// \brief Implementation of the RunActionMessenger class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunActionMessenger.hh"

#include "RunAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunActionMessenger::RunActionMessenger(RunAction* run)
    :G4UImessenger(),
     fRunAction(run),  fDetDir(0),
     fpixel_xCmd(0), fpixel_yCmd(0),fN_pixel_xCmd(0), fN_pixel_yCmd(0),
     fSource_typeCmd(0), fThresholdCmd(0), fFilenameCmd(0), fSaveProcPartCmd(0)
{
    fDetDir = new G4UIdirectory("/pixel/run/");
    fDetDir->SetGuidance("detector construction commands");

    fN_pixel_xCmd = new G4UIcmdWithAnInteger("/pixel/run/setN_pixel_x",this);
    fN_pixel_xCmd->SetGuidance("Set number of pixels to X.");
    fN_pixel_xCmd->SetParameterName("choice",false);
    fN_pixel_xCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fN_pixel_yCmd = new G4UIcmdWithAnInteger("/pixel/run/setN_pixel_y",this);
    fN_pixel_yCmd->SetGuidance("Set number of pixels to Y.");
    fN_pixel_yCmd->SetParameterName("choice",false);
    fN_pixel_yCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fpixel_xCmd = new G4UIcmdWithADoubleAndUnit("/pixel/run/setpixel_x", this);
    fpixel_xCmd->SetGuidance("Set pixel length X.");
    fpixel_xCmd->SetUnitCategory("Length");
    fpixel_xCmd->SetParameterName("choice",false);
    fpixel_xCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fpixel_yCmd = new G4UIcmdWithADoubleAndUnit("/pixel/run/setpixel_y", this);
    fpixel_yCmd->SetGuidance("Set pixel length Y.");
    fpixel_yCmd->SetUnitCategory("Length");
    fpixel_yCmd->SetParameterName("choice",false);
    fpixel_yCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fSource_typeCmd = new G4UIcmdWithAString("/pixel/run/setSource",this);
    fSource_typeCmd->SetGuidance("Select Source type.");
    fSource_typeCmd->SetParameterName("choice",false);
    fSource_typeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fThresholdCmd = new G4UIcmdWithADoubleAndUnit("/pixel/run/setThreshold",this);
    fThresholdCmd ->SetGuidance("Set Energy Threshold.");
    fThresholdCmd->SetDefaultValue((G4double)0.*eV);
    // fThresholdCmd->SetDefaultUnit("eV");
    fThresholdCmd ->SetUnitCategory("Energy");
    fThresholdCmd->SetParameterName("choice",false);
    fThresholdCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fFilenameCmd = new G4UIcmdWithAString("/pixel/run/setFilename",this);
    fFilenameCmd->SetGuidance("Select Filename.");
    fFilenameCmd->SetParameterName("choice",false);
    fFilenameCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fSaveProcPartCmd = new G4UIcmdWithABool("/pixel/run/setSaveProcPart",this);
    fSaveProcPartCmd->SetGuidance("Save File Process_Particles_Name.");
    fSaveProcPartCmd->SetParameterName("choice",false);
    fSaveProcPartCmd->SetDefaultValue(true);
    fSaveProcPartCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunActionMessenger::~RunActionMessenger()
{
    delete fN_pixel_xCmd;
    delete fN_pixel_yCmd;
    delete fpixel_xCmd;
    delete fpixel_yCmd;

    delete fSource_typeCmd;
    delete fThresholdCmd;
    delete fFilenameCmd;
    delete fSaveProcPartCmd;

    delete fDetDir;
    delete fRdecayDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunActionMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{

    if (command == fN_pixel_xCmd )
    {
        fRunAction->SetN_pixel_x(fN_pixel_xCmd->GetNewIntValue(newValue));
    }

    if (command == fN_pixel_yCmd )
    {
        fRunAction->SetN_pixel_y(fN_pixel_yCmd->GetNewIntValue(newValue));
    }
    if (command == fpixel_xCmd )
    {
        fRunAction->Setpixel_x(fpixel_xCmd->GetNewDoubleValue(newValue));
    }

    if (command == fpixel_yCmd )
    {
        fRunAction->Setpixel_y(fpixel_yCmd->GetNewDoubleValue(newValue));
    }


    if (command == fSource_typeCmd )
    {
        fRunAction->SetSource(newValue);
    }
    if (command == fThresholdCmd )
    {
        fRunAction->SetThreshold(fThresholdCmd->GetNewDoubleValue(newValue));
    }
    if (command == fFilenameCmd )
    {
        fRunAction->SetFilename(newValue);
    }

    if (command == fSaveProcPartCmd )
    {
        fRunAction->SetSaveProcePartname(fSaveProcPartCmd->GetNewBoolValue(newValue));
    }


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
