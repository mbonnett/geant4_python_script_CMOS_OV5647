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
/// \file DetectorMessenger.hh
/// \brief Definition of the DetectorMessenger class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorMessenger_h
#define DetectorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class DetectorConstruction;
class RunAction;

class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger;
class G4UIcmdWithABool;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorMessenger: public G4UImessenger
{
public:

    DetectorMessenger(DetectorConstruction* );
    ~DetectorMessenger();

    virtual void SetNewValue(G4UIcommand*, G4String);

private:

    DetectorConstruction*     fDetector;

    G4UIdirectory*              fRdecayDir;
    G4UIdirectory*              fDetDir;

    //G4UIcmdWithABool* fSource_typeCmd;

    G4UIcmdWithADoubleAndUnit*  fpixel_xCmd;
    G4UIcmdWithADoubleAndUnit*  fpixel_yCmd;
    G4UIcmdWithADoubleAndUnit*  fpixel_zCmd;
    G4UIcmdWithADoubleAndUnit*  fpixels_gapCmd;

    G4UIcmdWithAnInteger*       fN_pixel_xCmd;
    G4UIcmdWithAnInteger*       fN_pixel_yCmd;

    G4UIcmdWithADoubleAndUnit*  fDistance_SDCmd;
    G4UIcmdWithAString*         fSource_typeCmd;


    G4UIcmdWithABool*			fDetec_substraCmd;
    G4UIcmdWithADoubleAndUnit*  fSi_substrat_thickCmd;
    G4UIcmdWithABool*			fDetec_shieldCmd;
    G4UIcmdWithADoubleAndUnit*  fDetc_shield_thickCmd;

    G4UIcmdWithABool*			fDetec_SiO2_layers;
    G4UIcmdWithADoubleAndUnit*  fSiO2_layers_thinkCmd;
    G4UIcmdWithADoubleAndUnit*  fz_lensCmd;
    G4UIcmdWithADoubleAndUnit*  colr_fltrCmd;

    G4UIcmdWithADoubleAndUnit* fthick_shiel_down_Cmd;
    G4UIcmdWithADoubleAndUnit* fthick_shiel_up_Cmd;

    RunAction*               fRunAction;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
