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
/// \file RunAction.hh
/// \brief Definition of the RunAction class

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "G4Accumulable.hh"
#include "globals.hh"

class RunActionMessenger;
class DetectorConstruction;

/// Run action class

class RunAction : public G4UserRunAction
{
public:
    RunAction();
    virtual ~RunAction();

    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);

    void AddTrackLength(G4double length);

    void CountEvent()
    {
        fGoodEvents += 1;
    };

public:
    void SetN_pixel_x(G4int value);
    void SetN_pixel_y(G4int value);
    void Setpixel_x(G4double value);
    void Setpixel_y(G4double value);

    void SetSource(G4String value);
    void SetThreshold(G4double value);
    void SetFilename(G4String value);

    void SetSaveProcePartname(G4bool value);

public:
    G4int GetN_pixel_x();
    G4int GetN_pixel_y();
    G4double Getpixel_x();
    G4double Getpixel_y();

    G4String GetSource();
    G4double GetThreshold();
    G4String GetFilename();

    G4bool GetSaveProcePartname();

private:

    DetectorConstruction* fDetecConstruc;
    RunActionMessenger* fRunActionMessenger;

    G4int   fN_pixel_x;
    G4int   fN_pixel_y;
    G4double fpixel_x;
    G4double fpixel_y;

    G4String    fSource;
    G4double    fThreshold;
    G4String    fFilename;

    G4bool	     fSaveProcPart;

private:
    G4Accumulable<G4int>    fGoodEvents;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
