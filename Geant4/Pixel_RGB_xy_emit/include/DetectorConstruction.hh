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
/// \file DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;

class DetectorMessenger;

/// Detector construction class to define materials and geometry.
///
/// Crystals are positioned in Ring, with an appropriate rotation matrix.
/// Several copies of Ring are placed in the full detector.

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    DetectorConstruction();
    virtual ~DetectorConstruction();

public:
    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();

    void Setpixel_x(G4double value);
    void Setpixel_y(G4double value);
    void Setpixel_z(G4double value);
    void Setpixels_gap(G4double value);

    void SetN_pixel_x(G4int value);
    void SetN_pixel_y(G4int value);

    void SetDistance_SD(G4double value);
    void SetSource(G4String value);

    void SetDetec_substra(G4bool value);
    void SetSi_substrat_thick(G4double value);

    void SetDetec_shield(G4bool value);
    void SetDetc_shield_thick(G4double value);

    void SetDetec_SiO2_layers(G4bool value);
    void SetSiO2_layers_think(G4double value);
    void Setz_lens(G4double value);
    void Setz_colr_fltr(G4double value);

    void Set_thick_shiel_down(G4double value);
    void Set_thick_shiel_up(G4double value);

public:

    G4double Getpixel_x();
    G4double Getpixel_y();
    G4double Getpixel_z();
    G4double Getpixels_gap();

    G4int GetN_pixel_x();
    G4int GetN_pixel_y();

    G4double GetDistance_SD();
    G4String GetSource();

    G4bool GetDetec_substra();
    G4double GetSi_substrat_thick();

    G4bool GetDetec_shield();
    G4double GetDetc_shield_thick();

    G4bool GetDetec_SiO2_layers();
    G4double GetSiO2_layers_think();
    G4double Getz_lens();
    G4double Getz_colr_fltr();

    G4double Get_thick_shiel_down();
    G4double Get_thick_shiel_up();

private:

    DetectorMessenger* fDetectorMessenger;

    G4double    fpixel_x;
    G4double    fpixel_y;
    G4double    fpixel_z;
    G4double    fpixels_gap;

    G4int   fN_pixel_x;
    G4int   fN_pixel_y;

    G4double    fDistance_SD;
    G4String    fSource;

    G4bool		fDetec_substra;
    G4double    fSi_substrat_thick;

    G4bool		fDetec_shield;
    G4double    fDetc_shield_thick;

    G4bool		fDetec_SiO2_layers;
    G4double    fSiO2_layers_think;
    G4double    fz_lens;
    G4double    colr_fltr;

    G4double    fthick_shiel_down;
    G4double    fthick_shiel_up;

private:
    void DefineMaterials();

    G4bool  fCheckOverlaps;


};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
