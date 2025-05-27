//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                miki0             *
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
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Ellipsoid.hh"
#include "G4UnionSolid.hh"

#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSDoseDeposit.hh"
// G4PSCellCharge
#include "G4PSCellCharge.hh"

#include "G4VisAttributes.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"


#include "G4UnitsTable.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"

#include "G4PVReplica.hh"
#include "G4SDChargedFilter.hh"
#include "G4PSTrackLength.hh"
#include "G4PSPassageTrackLength.hh"
#include "G4Colour.hh"

#include "DetectorMessenger.hh"
#include "G4SDParticleFilter.hh"
#include "G4PSNofSecondary.hh"
#include "G4PSMinKinEAtGeneration.hh"
#include "G4PSNofStep.hh"

#include "EnergySD.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
    : G4VUserDetectorConstruction(),
      fDetectorMessenger(0),
      fCheckOverlaps(true)
{

    fDistance_SD = 0.0*mm;
//    fSource = "Am241";
    fSource = "gamma";

    fpixel_x = 1.4*um;
    fpixel_y = 1.4*um;
    fpixel_z = 2*um;
    fpixels_gap = 0.0*um;

    fN_pixel_x = 2592;
    fN_pixel_y = 1944;

    fDetec_substra = false;
    fSi_substrat_thick = 0*um;

    fDetec_shield = false;
    fDetc_shield_thick = 0*um;

    fDetec_SiO2_layers = true;
    fSiO2_layers_think = 0.225*um;

	  fz_lens = 0.5*um; // 0.735 um
    colr_fltr = 0.79*um; // 0.9um

    fthick_shiel_up = 2794*um; // 380*um;//110*um ;
    fthick_shiel_down = 381*um; // 380*um;//110*um ;

    DefineMaterials();

    fDetectorMessenger = new DetectorMessenger(this);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
    delete fDetectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  /*
    G4NistManager* man = G4NistManager::Instance();

    G4bool isotopes = false;

    G4Element*  O = man->FindOrBuildElement("O", isotopes);
    G4Element*  N = man->FindOrBuildElement("N", isotopes);
    G4Element* Si = man->FindOrBuildElement("Si", isotopes);
    G4Element* Lu = man->FindOrBuildElement("Lu", isotopes);

    G4Element* H = man->FindOrBuildElement("H", isotopes);
    G4Element* C = man->FindOrBuildElement("C", isotopes);
*/
//  G4Material* Strontium = new G4Material("Strontium", 38, 87.62*g/mole, *g/cm3,kStateUndefined, 273.15*kelvin, 1.0*atmosphere );
/*
    G4Material* LSO = new G4Material("Lu2SiO5", 7.4*g/cm3, 3);
    LSO->AddElement(Lu, 2);
    LSO->AddElement(Si, 1);
    LSO->AddElement(O, 5);

    // PDMS SiOC2H6
    // -------------
        G4Material* PDMS = new G4Material("PDMS", 1.34*g/cm3, 4);
        PDMS->AddElement(Si, 1 );
        PDMS->AddElement(O, 1 );
        PDMS->AddElement(C, 2 );
        PDMS->AddElement(H, 6 );

// Materials from Combination
//    G4Material* Air = new G4Material("Air",  1.205*mg/cm3, 2, kStateUndefined, 293.15*kelvin, 1.0*atmosphere );
//    Air->AddElement(N, 1);
//    Air->AddElement(O, 0);default_mat

    // Materials SiO2
//    G4Material* SiO2 = new G4Material("SiO2",  2.200*mg/cm3, 2);
//    SiO2->AddElement(Si, 1);
//    SiO2->AddElement(O, 2);

/////////////////////////////////////////////////////////////////////////////////////

    new G4Material("Galactic", 1., 1.008 * g / mole, 1.e-25 * g / cm3, kStateGas, 2.73 * kelvin, 3.e-18 * pascal);
    new G4Material("Gold", 79., 196.97 * g / mole, 19.32 * g / cm3);
    G4NistManager* nistManager = G4NistManager::Instance();
    G4Material* StainlessSteel = new G4Material("StainlessSteel", 8.06 * g / cm3, 6);
    StainlessSteel->AddElement(nistManager->FindOrBuildElement("C"), 0.001);
    StainlessSteel->AddElement(nistManager->FindOrBuildElement("Si"), 0.007);
    StainlessSteel->AddElement(nistManager->FindOrBuildElement("Cr"), 0.18);
    StainlessSteel->AddElement(nistManager->FindOrBuildElement("Mn"), 0.01);
    StainlessSteel->AddElement(nistManager->FindOrBuildElement("Fe"), 0.712);
    StainlessSteel->AddElement(nistManager->FindOrBuildElement("Ni"), 0.09);
/////////////////////////////////////////////////////////////////////////////////////
*/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{

//  G4Material* ll  = G4Material::GetMaterial("Lu");

    // detector Parameters
    if(fDetec_SiO2_layers==false) {
        fSiO2_layers_think = 0*um;
    }



    G4double detector_x = fN_pixel_x*fpixel_x;
    G4double detector_y = fN_pixel_y*fpixel_y;

    //  crystal
    G4double pixel_dX = fpixel_x - fpixels_gap,   pixel_dY = fpixel_y - fpixels_gap;

    //
    G4NistManager* nist = G4NistManager::Instance();
    G4Material* default_mat = nist->FindOrBuildMaterial("G4_AIR");
    //G4Material* cryst_mat   = nist->FindOrBuildMaterial("Lu2SiO5");
    G4Material* cmos_mat   = nist->FindOrBuildMaterial("G4_Si");
    //G4Material* aluminum   = nist->FindOrBuildMaterial("G4_Al");
    G4Material* Plexiglass = nist->FindOrBuildMaterial("G4_PLEXIGLASS");

    G4Material* Polyetilen = nist->FindOrBuildMaterial("G4_POLYETHYLENE");///cover
  //  G4Material* Polycarbonat = nist->FindOrBuildMaterial("G4_POLYCARBONATE");////lentes

    // auto scintillator = G4Material::GetMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
//"G4_POLYCARBONATE" ///lentes
//////////////"G4_POLYACRYLONITRILE" ///telas
//"G4_POLYETHYLENE"/////envasecover
///////////////"G4_POLYOXYMETHYLENE"//plastico duro
////////////////"G4_POLYPROPYLENE"
///////"G4_POLYSTYRENE" ////////envases de yougur

    // define Elements
    //
    G4Element* H  = new G4Element("Hydrogen", "H",  1,    1.01*g/mole);
    G4Element* C  = new G4Element("Carbon",   "C",  6,    12.01*g/mole);
    G4Element* N  = new G4Element("Nitrogen", "N",  7,    14.01*g/mole);
    G4Element* O  = new G4Element("Oxygen",   "O",  8,    16.00*g/mole);
    G4Element* Si = new G4Element("Silicon",  "Si", 14.,  28.09*g/mole);
    //
    // define simple materials
    //G4Element* Si  =new G4Material("Silicon"  symbol="Si", z=14, a=28.09*g/mole, density= 2.330*g/cm3)
    /*
    G4Element* H  = manager->FindOrBuildElement(1);
    G4Element* C  = manager->FindOrBuildElement(6);
    G4Element* N  = manager->FindOrBuildElement(7);
    G4Element* O  = manager->FindOrBuildElement(8);
    G4Element* Si = manager->FindOrBuildElement(14);
    */
  // PDMS SiOC2H6
  // -------------
    G4Material* PDMS = new G4Material("PDMS", 1.34*g/cm3, 4);
    PDMS->AddElement(Si, 1 );
    PDMS->AddElement(O, 1 );
    PDMS->AddElement(C, 2 );
    PDMS->AddElement(H, 6 );

    ///Photoresist (C10H6N2O)
    G4Material* PhotoR = new G4Material("PhotoR", 1.29*g/cm3, 4);
    PhotoR->AddElement(C, 10 );
    PhotoR->AddElement(H, 6 );
    PhotoR->AddElement(N, 2 );
    PhotoR->AddElement(O, 1 );

    // PMMA C5H8O2 ( Acrylic )
    // -------------
    G4Material* Acrylic = new G4Material("Acrylic", 1.19*g/cm3, 3);
    Acrylic->AddElement(C,  5);
    Acrylic->AddElement(H, 8);
    Acrylic->AddElement(O, 2);


    //Epoxy (for FR4 )
    G4Material* Epoxy = new G4Material("Epoxy" , 1.2*g/cm3, 2);
    Epoxy->AddElement(H, 2);
    Epoxy->AddElement(C, 2);


    //G4Material* Brass = nist->FindOrBuildMaterial("G4_BRASS");
    G4Material* Stainless_steel = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
    G4Material* copper = nist->FindOrBuildMaterial("G4_Cu");
    G4Material* Beryllium = nist->FindOrBuildMaterial("G4_Be");

    G4Material* SiO2 = nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");

    G4Material* Cellophane = nist->FindOrBuildMaterial("G4_CELLULOSE_CELLOPHANE");
    // Material: Vacuum
  //  G4Material* Vacuum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
    //G4Material* default_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");

  //  if(fSource=="e" ){ default_mat=Vacuum;}
  //  else{ default_mat=Air;}

      // Material: Vacuum
    //  G4Material* Vacuum = new G4Material("Vacuum", 1.0 , 1.01*g/mole, 1.0E-25*g/cm3, kStateGas, 2.73*kelvin, 3.0E-18*pascal );

    //Source shielding_ Sr90 Cs137
    G4double thick_z =  380*um;//110*um ;
    G4double ring_R2 = 12700 *um, ring_R1 = 3175*um;
    G4double Source_Height = 3175*um - fthick_shiel_up - fthick_shiel_down;
    G4double epoxy_Height = fthick_shiel_up;
    G4double thick_paper =  160*um;//110*um ;

    //Source shielding_ Fe55
    G4double monel_thick_z = 250*um, be_thick_z = 250*um;
    G4double monel_z = 2.5*mm - monel_thick_z;
//Diámetro total: 1 pulgada (25.4 mm)
//Diámetro activo: 0.197 pulgada (5 mm)  or 0.25in  6350*um
//Altura: 0.125 pulgada (3.18 mm)   3175*um

    G4double place_Detec_Z;

    //
    // World
    //
    G4double world_sizeX = 5*cm;
    G4double world_sizeY = 5*cm;
    G4double world_sizeZ = 10*cm;


// Visualization attributes


    G4VisAttributes * tomato= new G4VisAttributes( G4Colour(255/255.,99/255.,71/255. ));
    tomato -> SetVisibility(true);
    G4VisAttributes * silver= new G4VisAttributes( G4Colour(192/255.,192/255.,192/255. ));
    silver -> SetVisibility(true);
    G4VisAttributes * goldenrod= new G4VisAttributes( G4Colour(218/255.,165/255.,32/255. ));
    goldenrod -> SetVisibility(true);
    G4VisAttributes * orange= new G4VisAttributes( G4Colour(255/255.,165/255.,0/255. ));
    orange -> SetVisibility(true);
    G4VisAttributes * green= new G4VisAttributes( G4Colour(0/255.,128/255.,0/255. ));
    green -> SetVisibility(true);
    G4VisAttributes * yellow= new G4VisAttributes( G4Colour(255/255., 255/255. ,0/255. ));
    yellow -> SetVisibility(true);
    G4VisAttributes * brown= new G4VisAttributes( G4Colour(153/255.,102/255.,0/255. ));
    brown -> SetVisibility(true);

    G4VisAttributes * blue= new G4VisAttributes( G4Colour(0/255., 0/255., 255/255. ));
    blue -> SetVisibility(true);




    G4Box* solidWorld =
        new G4Box("World",                       //its name
                  0.5*world_sizeX, 0.5*world_sizeY, 0.5*world_sizeZ); //its size

    G4LogicalVolume* logicWorld =
        new G4LogicalVolume(solidWorld,          //its solid
                            default_mat,         //its material
                            "World");            //its name

    G4VPhysicalVolume* physWorld =
        new G4PVPlacement(0,                     //no rotation
                          G4ThreeVector(),       //at (0,0,0)
                          logicWorld,            //its logical volume
                          "World",               //its name
                          0,                     //its mother  volume
                          false,                 //no boolean operation
                          0,                     //copy number
                          fCheckOverlaps);       // checking overlaps

/////////////////////////////////

if(fSource=="e" ||fSource=="alpha")
{
    place_Detec_Z=fDistance_SD+fpixel_z/2+fSiO2_layers_think;
}


    if(fSource=="Sr90" ||fSource=="Cs137" )
    {

        place_Detec_Z=fDistance_SD+Source_Height/2+fpixel_z/2+fthick_shiel_down+fz_lens+colr_fltr+fSiO2_layers_think;

        //Source shielding
//
        G4Tubs* solid_Sshield_out = new G4Tubs("Sshield_out", //name
                                               ring_R1, ring_R2, Source_Height/2, 0.,twopi); //dimensions

        G4LogicalVolume* logic_Sshield_out =
            new G4LogicalVolume(solid_Sshield_out,           //its solid
                                Plexiglass,         //its material
                                "Sshield_out");             //its name


        new G4PVPlacement(0,                         //no rotation
                          G4ThreeVector(0,0,0),             //at (0,0,0)
                          logic_Sshield_out,                //logical volume
                          "Sshield_out",                    //name
                          logicWorld,                      //mother  volume
                          false,                       //no boolean operation
                          0,                       //copy number
                          fCheckOverlaps);         // checking overlaps
        logic_Sshield_out->SetVisAttributes(tomato);

        //cap Epoxy shield
        G4Tubs* solid_epoxy_Sshield_down = new G4Tubs("epoxy_Sshield_down", //name
                                               ring_R1/2, ring_R1, Source_Height/2, 0.,twopi); //dimensions

        G4LogicalVolume* logic_epoxy_Sshield_down =
            new G4LogicalVolume(solid_epoxy_Sshield_down,           //its solid
                                Epoxy,         //its material
                                "epoxy_Sshield_down");             //its name


        new G4PVPlacement(0,                         //no rotation
                          G4ThreeVector(0,0,0),             //at (0,0,0)
                          logic_epoxy_Sshield_down,                //logical volume
                          "epoxy_Sshield_down",                    //name
                          logicWorld,                      //mother  volume
                          false,                       //no boolean operation
                          0,                       //copy number
                          fCheckOverlaps);         // checking overlaps
        logic_epoxy_Sshield_down->SetVisAttributes(silver);

        G4Tubs* solid_epoxy_Sshield = new G4Tubs("epoxy_Sshield", //name
                0, ring_R1, epoxy_Height/2, 0.,twopi); //dimensions

        G4LogicalVolume* logic_epoxy_Sshield =
            new G4LogicalVolume(solid_epoxy_Sshield,           //its solid
                                Epoxy,         //its material
                                "epoxy_Sshield");             //its name


        new G4PVPlacement(0,                         //no rotation
                          G4ThreeVector(0,0,-epoxy_Height/2-Source_Height/2),             //at (0,0,0)
                          logic_epoxy_Sshield,                //logical volume
                          "epoxy_Sshield",                    //name
                          logicWorld,                      //mother  volume
                          false,                       //no boolean operation
                          0,                       //copy number
                          fCheckOverlaps);         // checking overlaps
        logic_epoxy_Sshield->SetVisAttributes(silver);

        //cap Source tub shieldingup
       if(fthick_shiel_up>0*um)
        {
        G4Tubs* solid_tub_Sshield_up = new G4Tubs("tub_Sshield_up", //name
                ring_R1, ring_R2, fthick_shiel_up/2, 0.,twopi); //dimensions

        G4LogicalVolume* logic_tub_Sshield_up =
            new G4LogicalVolume(solid_tub_Sshield_up,           //its solid
                                Plexiglass,         //its material
                                "tub_Sshield_up");             //its name


        new G4PVPlacement(0,                         //no rotation
                          G4ThreeVector(0,0,-Source_Height/2-fthick_shiel_up/2),             //at (0,0,0)
                          logic_tub_Sshield_up,                //logical volume
                          "tub_Sshield_up",                    //name
                          logicWorld,                      //mother  volume
                          false,                       //no boolean operation
                          0,                       //copy number
                          fCheckOverlaps);         // checking overlaps
        logic_tub_Sshield_up->SetVisAttributes(tomato);
        }

//cap paper shieldingup
        G4Tubs* solid_paper_Sshield_up = new G4Tubs("paper_Sshield_up", //name
                0, ring_R2, thick_paper/2, 0.,twopi); //dimensions

        G4LogicalVolume* logic_paper_Sshield_up =
            new G4LogicalVolume(solid_paper_Sshield_up,           //its solid
                                Cellophane,         //its material
                                "paper_Sshield_up");             //its name


        new G4PVPlacement(0,                         //no rotation
                          G4ThreeVector(0,0,-fthick_shiel_up-Source_Height/2-thick_paper/2),             //at (0,0,0)
                          logic_paper_Sshield_up,                //logical volume
                          "paper_Sshield_up",                    //name
                          logicWorld,                      //mother  volume
                          false,                       //no boolean operation
                          0,                       //copy number
                          fCheckOverlaps);         // checking overlaps
        logic_paper_Sshield_up->SetVisAttributes(yellow);


        //cap Source shielding down
       if(fthick_shiel_down>0*um)
       {
        G4Tubs* solid_cap_Sshield_down = new G4Tubs("cap_Sshield_down", //name
                0, ring_R2, fthick_shiel_down/2, 0.,twopi); //dimensions

        G4LogicalVolume* logic_cap_Sshield_down =
            new G4LogicalVolume(solid_cap_Sshield_down,           //its solid
                                Plexiglass,         //its material
                                "cap_Sshield_down");             //its name
        logic_cap_Sshield_down->SetVisAttributes(tomato);

        new G4PVPlacement(0,                         //no rotation
                          G4ThreeVector(0,0,1*Source_Height/2+fthick_shiel_down/2),             //at (0,0,0)
                          logic_cap_Sshield_down,                //logical volume
                          "cap_Sshield_down",                    //name
                          logicWorld,                      //mother  volume
                          false,                       //no boolean operationGetNewBoolValue
                          0,                       //copy number
                          fCheckOverlaps);         // checking overlaps
        }

    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    if(fSource=="Am241" )
    {
        //Source shielding_ Am
        ring_R1= 3000*um, Source_Height=2000*um;
        G4double Source_Height_0=500*um, ring_R0= 2000*um;

        place_Detec_Z=fDistance_SD+Source_Height/2+Source_Height_0/2+fpixel_z/2+fz_lens+colr_fltr+fSiO2_layers_think;

        //Source shielding_ Am241


        //Source shielding
//
        G4Tubs* solid_Sshield_out = new G4Tubs("Sshield_out", //name
                                               ring_R1, ring_R2, (Source_Height-Source_Height_0)/2, 0.,twopi); //dimensions

        G4LogicalVolume* logic_Sshield_out =
            new G4LogicalVolume(solid_Sshield_out,           //its solid
                                Polyetilen,         //its material
                                "Sshield_out");             //its name


        new G4PVPlacement(0,                         //no rotation
                          G4ThreeVector(0,0,Source_Height_0),             //at (0,0,0)
                          logic_Sshield_out,                //logical volume
                          "Sshield_out",                    //name
                          logicWorld,                      //mother  volume
                          false,                       //no boolean operation
                          0,                       //copy number
                          fCheckOverlaps);         // checking overlaps
        logic_Sshield_out->SetVisAttributes(tomato);
logic_Sshield_out->SetVisAttributes (G4VisAttributes::GetInvisible());

        //cap Source shielding_up

        G4Tubs* solid_cap_Sshield_up = new G4Tubs("cap_Sshield_up", //name
                Source_Height, ring_R2, Source_Height_0/2, 0.,twopi); //dimensions

        G4LogicalVolume* logic_cap_Sshield_up =
            new G4LogicalVolume(solid_cap_Sshield_up,           //its solid
                                Polyetilen,         //its material
                                "cap_Sshield_up");             //its name


        new G4PVPlacement(0,                         //no rotation
                          G4ThreeVector(0,0,-Source_Height/2+Source_Height_0),             //at (0,0,0)
                          logic_cap_Sshield_up,                //logical volume
                          "cap_Sshield_up",                    //name
                          logicWorld,                      //mother  volume
                          false,                       //no boolean operation
                          0,                       //copy number
                          fCheckOverlaps);         // checking overlaps
        logic_cap_Sshield_up->SetVisAttributes(tomato);
logic_cap_Sshield_up->SetVisAttributes (G4VisAttributes::GetInvisible());

        //cap Source shielding_top_up

        G4Tubs* solid_cap_Sshield_top_up = new G4Tubs("cap_Sshield_up", //name
                0, ring_R2, thick_z/2, 0.,twopi); //dimensions

        G4LogicalVolume* logic_cap_Sshield_top_up =
            new G4LogicalVolume(solid_cap_Sshield_top_up,           //its solid
                                Polyetilen,         //its material
                                "cap_Sshield_top_up");             //its name


        new G4PVPlacement(0,                         //no rotation
                          G4ThreeVector(0,0,-Source_Height/2+Source_Height_0/2-thick_z/2),             //at (0,0,0)
                          logic_cap_Sshield_top_up,                //logical volume
                          "cap_Sshield_top_up",                    //name
                          logicWorld,                      //mother  volume
                          false,                       //no boolean operation
                          0,                       //copy number
                          fCheckOverlaps);         // checking overlaps
        logic_cap_Sshield_top_up->SetVisAttributes(tomato);
//logic_cap_Sshield_top_up->SetVisAttributes (G4VisAttributes::GetInvisible());

        /*
                //cap Source shielding down

                G4Tubs* solid_cap_Sshield_down = new G4Tubs("cap_Sshield_down", //name
                        0, ring_R2, thick_z/2, 0.,twopi); //dimensions

                G4LogicalVolume* logic_cap_Sshield_down =
                    new G4LogicalVolume(solid_cap_Sshield_down,           //its solid
                                        Polyetilen,         //its material
                                        "cap_Sshield_down");             //its name
                logic_cap_Sshield_down->SetVisAttributes(tomato);

                new G4PVPlacement(0,                         //no rotation
                                  G4ThreeVector(0,0,1*Source_Height/2+thick_z/2),             //at (0,0,0)
                                  logic_cap_Sshield_down,                //logical volume
                                  "cap_Sshield_down",                    //name
                                  logicWorld,                      //mother  volume
                                  false,                       //no boolean operationGetNewBoolValue
                                  0,                       //copy number
                                  fCheckOverlaps);         // checking overlaps
        */
        //cap Source Stainless_steel

        G4Tubs* solid_cap_Stainless_steel = new G4Tubs("cap_Sshield_up", //name
                ring_R0, ring_R1, (Source_Height-Source_Height_0)/2, 0.,twopi); //dimensions

        G4LogicalVolume* logic_cap_Stainless_steel =
            new G4LogicalVolume(solid_cap_Stainless_steel,           //its solid
                                Stainless_steel,         //its material
                                "cap_Stainless_steel");             //its name


        new G4PVPlacement(0,                         //no rotation
                          G4ThreeVector(0,0,Source_Height_0),             //at (0,0,0)
                          logic_cap_Stainless_steel,                //logical volume
                          "cap_Stainless_steel",                    //name
                          logicWorld,                      //mother  volume
                          false,                       //no boolean operation
                          0,                       //copy number
                          fCheckOverlaps);         // checking overlaps
        logic_cap_Stainless_steel->SetVisAttributes(silver);



    }

/////////////////////////////////////////////////////////////////////////////////////////////////////////


    if(fSource=="Fe55")
    {

        place_Detec_Z=fDistance_SD+be_thick_z+monel_thick_z/2+fpixel_z/2+fz_lens+colr_fltr+fSiO2_layers_think; //for Fe55

        //Source shielding_ Fe55
//


        G4Tubs* solid_Sshield_out = new G4Tubs("Sshield_out", //name
                                               6*mm, 7.5*mm, monel_z/1, 0.,twopi); //dimensions

        G4LogicalVolume* logic_Sshield_out =
            new G4LogicalVolume(solid_Sshield_out,           //its solid
                                Stainless_steel,         //its material
                                "Sshield_out");             //its name


        new G4PVPlacement(0,                         //no rotation
                          G4ThreeVector(0,0,-1*monel_z+1.5*monel_thick_z),             //at (0,0,0)
                          logic_Sshield_out,                //logical volume
                          "Sshield_out",                    //name
                          logicWorld,                      //mother  volume
                          false,                       //no boolean operation
                          0,                       //copy number
                          fCheckOverlaps);         // checking overlaps

        //cap Source shieldingup

        G4Tubs* solid_cap_Sshield_up = new G4Tubs("cap_Sshield_up", //name
                0, 7.5*mm, monel_thick_z/2, 0.,twopi); //dimensions

        G4LogicalVolume* logic_cap_Sshield_up =
            new G4LogicalVolume(solid_cap_Sshield_up,           //its solid
                                Stainless_steel,         //its material
                                "cap_Sshield_up");             //its name


        new G4PVPlacement(0,                         //no rotation
                          G4ThreeVector(0,0,-2*monel_z+monel_thick_z),             //at (0,0,0)
                          logic_cap_Sshield_up,                //logical volume
                          "cap_Sshield_up",                    //name
                          logicWorld,                      //mother  volume
                          false,                       //no boolean operation
                          0,                       //copy number
                          fCheckOverlaps);         // checking overlaps

        //cap Source shielding down

        G4Tubs* solid_cap_Sshield_down = new G4Tubs("cap_Sshield_down", //name
                5*mm, 7.5*mm, monel_thick_z/2, 0.,twopi); //dimensions

        G4LogicalVolume* logic_cap_Sshield_down =
            new G4LogicalVolume(solid_cap_Sshield_down,           //its solid
                                Stainless_steel,         //its material
                                "cap_Sshield_down");             //its name


        new G4PVPlacement(0,                         //no rotation
                          G4ThreeVector(0,0,be_thick_z+monel_thick_z),             //at (0,0,0)
                          logic_cap_Sshield_down,                //logical volume
                          "cap_Sshield_down",                    //name
                          logicWorld,                      //mother  volume
                          false,                       //no boolean operation
                          0,                       //copy number
                          fCheckOverlaps);         // checking overlaps

//cap be win

        G4Tubs* solid_be_win = new G4Tubs("be_win", //name
                                          0, 6*mm, be_thick_z/2, 0.,twopi); //dimensions

        G4LogicalVolume* logic_be_win =
            new G4LogicalVolume(solid_be_win,           //its solid
                                Beryllium,         //its material
                                "be_win");             //its name

        logic_be_win->SetVisAttributes(silver);

        new G4PVPlacement(0,                         //no rotation
                          G4ThreeVector(0,0,be_thick_z),             //at (0,0,0)
                          logic_be_win,                //logical volume
                          "be_win",                    //name
                          logicWorld,                      //mother  volume
                          false,                       //no boolean operation
                          0,                       //copy number
                          fCheckOverlaps);         // checking overlaps


//cap copper_shield


        G4Tubs* solid_cu_shield = new G4Tubs("cu_shield", //name
                                             0, 6*mm, monel_z-0.5*be_thick_z-monel_thick_z/2, 0.,twopi); //dimensions

        G4LogicalVolume* logic_cu_shield =
            new G4LogicalVolume(solid_cu_shield,           //its solid
                                copper,         //its material
                                "cu_shield");             //its name

        logic_cu_shield->SetVisAttributes(goldenrod);


        new G4PVPlacement(0,                         //no rotation
                          G4ThreeVector(0,0,-1*monel_z+0.5*be_thick_z),             //at (0,0,0)
                          logic_cu_shield,                //logical volume
                          "cu_shield",                    //name
                          logicWorld,                      //mother  volume
                          false,                       //no boolean operation
                          0,                       //copy number
                          fCheckOverlaps);         // checking overlaps

    }

////////////////////////////////////////////////////////////////////////////////

/////////Lens---PDMS
    if(fz_lens>0*um )
    {
    G4Box* solidDetec_lens =
        new G4Box("Detec_lens",                       //its name
                  0.5*detector_x, 0.5*detector_y, 0.5*fz_lens); //its size

    G4LogicalVolume* logicDetec_lens =
        new G4LogicalVolume(solidDetec_lens,       //its solid
                            default_mat,         //its material
                            "Detec_lens");         //its name

    //  G4VSolid*
  G4ThreeVector zTrans(0, 0, 0.25*fz_lens);
//  auto* solidLen1 = new G4Ellipsoid("lent1",0.5*pixel_dX, 0.5*pixel_dY, 0.5*fz_lens, -0.5*fz_lens, 0);
//  auto* solidLen2 = new G4Box("lent2",0.5*pixel_dX, 0.5*pixel_dY, 0.25*fz_lens);
  auto* solidLen1 = new G4Ellipsoid("lent1",0.5*fpixel_x, 0.5*fpixel_y, 0.5*fz_lens, -0.5*fz_lens, 0);
  auto* solidLen2 = new G4Box("lent2",0.5*fpixel_x, 0.5*fpixel_y, 0.25*fz_lens);
  G4UnionSolid* solidLen =  new G4UnionSolid("lent", solidLen1, solidLen2, 0, zTrans);


  G4LogicalVolume* logicLen =
              new G4LogicalVolume(solidLen,          //its solid
                                  PDMS,           //its material
                                  "LenLV");        //its name

logicLen->SetVisAttributes(tomato);
///
// place len
//
    for (G4int ipxl_x= 0; ipxl_x < fN_pixel_x ; ipxl_x++)
    {
        for (G4int ipxl_y= 0; ipxl_y < fN_pixel_y ; ipxl_y++)
        {

            G4double position_x = ((2*ipxl_x+1)*fpixel_x-detector_x)/2;
            G4double position_y = ((2*ipxl_y+1)*fpixel_y-detector_y)/2;

            new G4PVPlacement(0,                       //no rotation
                              G4ThreeVector(position_x,position_y,0.),//at (0,0,0)
                              logicLen,            //its logical volume
                              "Detec_lens",             //its name
                              logicDetec_lens,             //its mother  volume
                              false,                 //no boolean operation
                              ipxl_x + ipxl_y*fN_pixel_x+1,                 //copy number
                              0);       // checking overlaps
        }
    }

    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(0,0,1.*place_Detec_Z-fz_lens/2-colr_fltr-fSiO2_layers_think-0.5*fpixel_z),       //at (0,0,0)
                      logicDetec_lens,            //its logical volume
                      "Detec_lens",               //its name
                      logicWorld,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      fCheckOverlaps);       // checking overlaps

    //logicDetec_lens->SetVisAttributes(tomato);
    logicDetec_lens->SetVisAttributes (G4VisAttributes::GetInvisible());
    }
////////////////////////////////////////////////////////////////////////////////

/////////col0r_filter---Photoresist (C10H6N2O)
    if(colr_fltr>0*um )
    {
    G4Box* solidDetec_colr_fltr =
        new G4Box("Detec_colr_fltr",                       //its name
                  0.5*detector_x, 0.5*detector_y, 0.5*colr_fltr); //its size

    G4LogicalVolume* logicDetec_colr_fltr =
        new G4LogicalVolume(solidDetec_colr_fltr,       //its solid
                            PhotoR,         //its material
                            "Detec_colr_fltr");         //its name

// define filter

//G4Box* solidFiltr = new G4Box("filter", 0.5*pixel_dX, 0.5*pixel_dY, 0.5*colr_fltr);
G4Box* solidFiltr = new G4Box("filter", 0.5*fpixel_x, 0.5*fpixel_y, 0.5*colr_fltr);

G4LogicalVolume* logicFiltr =
    new G4LogicalVolume(solidFiltr,          //its solid
                        PhotoR,           //its material
                        "FilterLV");        //its name

logicFiltr->SetVisAttributes(yellow);

// place filter
//
for (G4int ipxl_x= 0; ipxl_x < fN_pixel_x ; ipxl_x++)
{
    for (G4int ipxl_y= 0; ipxl_y < fN_pixel_y ; ipxl_y++)
    {

        G4double position_x = ((2*ipxl_x+1)*fpixel_x-detector_x)/2;
        G4double position_y = ((2*ipxl_y+1)*fpixel_y-detector_y)/2;

        new G4PVPlacement(0,                       //no rotation
                          G4ThreeVector(position_x,position_y,0.),//at (0,0,0)
                          logicFiltr,            //its logical volume
                          "Detec_colr_fltr",             //its name
                          logicDetec_colr_fltr,             //its mother  volume
                          false,                 //no boolean operation
                          ipxl_x + ipxl_y*fN_pixel_x+1,                 //copy number
                          0);       // checking overlaps
    }
}

    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(0, 0, 1.*place_Detec_Z-0.5*colr_fltr-fSiO2_layers_think-0.5*fpixel_z),       //at (0,0,0)
                      logicDetec_colr_fltr,            //its logical volume
                      "Detec_colr_fltr",               //its name
                      logicWorld,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      fCheckOverlaps);       // checking overlaps

    logicDetec_colr_fltr->SetVisAttributes(yellow);
    }
////////////////////////////////////////////////////////////////////////////////


/////////SiO2///////////////////////////////////////////////////////////////////
    if(fDetec_SiO2_layers==true)
    {
        G4Box* solidDetec_SiO2_layers =
            new G4Box("Detec_SiO2_layers",                       //its name
                      0.5*detector_x, 0.5*detector_y, 0.5*fSiO2_layers_think); //its size

        G4LogicalVolume* logicDetec_SiO2_layers =
            new G4LogicalVolume(solidDetec_SiO2_layers,       //its solid
                                SiO2,         //its material
                                "Detec_SiO2_layers");         //its name

        new G4PVPlacement(0,                     //no rotation
                          G4ThreeVector(0,0,1.*place_Detec_Z-fSiO2_layers_think/2-0.5*fpixel_z),       //at (0,0,0)
                          logicDetec_SiO2_layers,            //its logical volume
                          "Detec_SiO2_layers",               //its name
                          logicWorld,                     //its mother  volume
                          false,                 //no boolean operation
                          0,                     //copy number
                          fCheckOverlaps);       // checking overlaps

        logicDetec_SiO2_layers->SetVisAttributes(blue);
    }
////////////////////////////////////////////////////////////////////////////////

////Detector

    G4Box* solidDetector =
        new G4Box("Detector",                       //its name
                  0.5*detector_x, 0.5*detector_y, 0.5*fpixel_z); //its size

    G4LogicalVolume* logicDetector =
        new G4LogicalVolume(solidDetector,       //its solid
                            cmos_mat,         //its material
                            "Detector");         //its name
    // define crystal

    G4Box* solidCryst = new G4Box("crystal", 0.5*pixel_dX, 0.5*pixel_dY, 0.5*fpixel_z);

    G4LogicalVolume* logicCryst =
        new G4LogicalVolume(solidCryst,          //its solid
                            cmos_mat,           //its material
                            "CrystalLV");        //its name

    logicCryst->SetVisAttributes(green);

    // place Detector
    //
    for (G4int ipxl_x= 0; ipxl_x < fN_pixel_x ; ipxl_x++)
    {
        for (G4int ipxl_y= 0; ipxl_y < fN_pixel_y ; ipxl_y++)
        {

            G4double position_x = ((2*ipxl_x+1)*fpixel_x-detector_x)/2;
            G4double position_y = ((2*ipxl_y+1)*fpixel_y-detector_y)/2;

            new G4PVPlacement(0,                       //no rotation
                              G4ThreeVector(position_x,position_y,0.),//at (0,0,0)
                              logicCryst,            //its logical volume
                              "Detector",             //its name
                              logicDetector,             //its mother  volume
                              false,                 //no boolean operation
                              ipxl_x + ipxl_y*fN_pixel_x+1,                 //copy number
                              0);       // checking overlaps
        }
    }

    new G4PVPlacement(0,                       //no rotation
                      G4ThreeVector(0.,0.,place_Detec_Z), //at (0,0,0)
                      logicDetector,           //its logical volume
                      "Detector",              //its name
                      logicWorld,              //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      fCheckOverlaps);         // checking overlaps

    //logicDetector->SetVisAttributes(green);
    logicDetector->SetVisAttributes (G4VisAttributes::GetInvisible());

/////////SUBSTRAT///////////////////////////////////////////////////////////////
    if(fDetec_substra==true)
    {
        G4Box* solidDetecBase_Si =
            new G4Box("DetecBase_Si",                       //its name
                      0.5*detector_x, 0.5*detector_y, 0.5*fSi_substrat_thick); //its size

        G4LogicalVolume* logicDetecBase_Si =
            new G4LogicalVolume(solidDetecBase_Si,       //its solid
                                cmos_mat,         //its material
                                "DetecBase_Si");         //its name

        new G4PVPlacement(0,                     //no rotation
                          G4ThreeVector(0,0,1.*place_Detec_Z+0.5*fSi_substrat_thick+0.5*fpixel_z),       //at (0,0,0)
                          logicDetecBase_Si,            //its logical volume
                          "DetecBase_Si",               //its name
                          logicWorld,                     //its mother  volume
                          false,                 //no boolean operation
                          0,                     //copy number
                          fCheckOverlaps);       // checking overlaps

        logicDetecBase_Si->SetVisAttributes(brown);
    }

    if(fDetec_shield==true)
    {
        G4Box* solidDetecBase_Ss =
            new G4Box("DetecBase_Ss",                       //its name
                      0.5*detector_x, 0.5*detector_y, 0.5*fDetc_shield_thick); //its size

        G4LogicalVolume* logicDetecBase_Ss =
            new G4LogicalVolume(solidDetecBase_Ss,       //its solid
                                Plexiglass,         //its material
                                "DetecBase_Ss");         //its name

        new G4PVPlacement(0,                     //no rotation
                          G4ThreeVector(0,0,1.*place_Detec_Z+0.5*fDetc_shield_thick+1*fSi_substrat_thick+0.5*fpixel_z),       //at (0,0,0)
                          logicDetecBase_Ss,            //its logical volume
                          "DetecBase_Ss",               //its name
                          logicWorld,                     //its mother  volume
                          false,                 //no boolean operation
                          0,                     //copy number
                          fCheckOverlaps);       // checking overlaps

        //logicDetecBase_Ss->SetVisAttributes(goldenrod);



        G4Box* solidDetecBase0 =
            new G4Box("DetecBase0",                       //its name
                      0.55*detector_x, 0.55*detector_y, 0.5*(fpixel_z+fSi_substrat_thick+fDetc_shield_thick)); //its size

        G4Box* solidDetecBase1 =
            new G4Box("DetecBase1",                       //its name
                      0.5*detector_x, 0.5*detector_y, 0.75*fpixel_z+0.5*fSi_substrat_thick+0.5*fDetc_shield_thick); //its size

        G4SubtractionSolid* solidDetecBase2 =
            new G4SubtractionSolid("DetecBase", solidDetecBase0, solidDetecBase1);

        G4LogicalVolume* logicDetecBase2 =
            new G4LogicalVolume(solidDetecBase2,       //its solid
                                Plexiglass,         //its material
                                "DetecBase");         //its name

        new G4PVPlacement(0,                     //no rotation
                          G4ThreeVector(0,0,place_Detec_Z+0*fpixel_z+0.5*(fSi_substrat_thick+fDetc_shield_thick)),       //at (0,0,0)
                          logicDetecBase2,            //its logical volume
                          "DetecBase",               //its name
                          logicWorld,                     //its mother  volume
                          false,                 //no boolean operation
                          0,                     //copy number
                          fCheckOverlaps);       // checking overlaps

    }

///////////////////////////////////////////////////////////////////////////////////////

    // logicDetecBase->SetVisAttributes(tomato);

//
    // Visualization attributes
    //
    // logicDetector->SetVisAttributes (G4VisAttributes::GetInvisible());
    // Detector_pixel_LV->SetVisAttributes (G4VisAttributes::GetInvisible());

    logicWorld->SetVisAttributes (G4VisAttributes::GetInvisible());

    auto simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
    simpleBoxVisAtt->SetVisibility(true);
    logicDetector->SetVisAttributes(simpleBoxVisAtt);

    // Print materials
    /*    G4cout << *(G4Material::GetMaterialTable()) << G4endl;
    //////////////////////////////////////////////////////////////////////////
        // print parameters
        //
        G4cout
                << G4endl
                << "------------------------------------------------------------" << G4endl
                << "---> The Detector is " << fN_pixel_x*fN_pixel_y << " pixel of: [ "
                << fpixel_x/um << "um of " << cmos_mat->GetName()
                << " with a gap "
                << fpixels_gap/um << "um "<< " ] " << G4endl
                << "------------------------------------------------------------" << G4endl;
    */
    //
    //always return the physical World
    //
    return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::Setpixel_x(G4double value)
{
    fpixel_x = value;
    //G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double DetectorConstruction::Getpixel_x()
{
    return fpixel_x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::Setpixel_y(G4double value)
{
    fpixel_y = value;
    //G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double DetectorConstruction::Getpixel_y()
{
    return fpixel_y;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::Setpixel_z(G4double value)
{
    fpixel_z = value;
    //G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double DetectorConstruction::Getpixel_z()
{
    return fpixel_z;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::Setpixels_gap(G4double value)
{
    fpixels_gap = value;
    //G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double DetectorConstruction::Getpixels_gap()
{
    return fpixels_gap;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetN_pixel_x(G4int value)
{
    fN_pixel_x = value;
    //G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4int DetectorConstruction::GetN_pixel_x()
{
    return fN_pixel_x;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetN_pixel_y(G4int value)
{
    fN_pixel_y = value;
    //G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4int DetectorConstruction::GetN_pixel_y()
{
    return fN_pixel_y;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetDistance_SD(G4double value)
{
    fDistance_SD = value;
    //G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double DetectorConstruction::GetDistance_SD()
{
    return fDistance_SD;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetSource(G4String value)
{
    fSource = value;
    //G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4String DetectorConstruction::GetSource()
{
    return fSource;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetDetec_substra(G4bool value)
{
    fDetec_substra = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool DetectorConstruction::GetDetec_substra()
{
    return fDetec_substra;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetSi_substrat_thick(G4double value)
{
    fSi_substrat_thick = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double DetectorConstruction::GetSi_substrat_thick()
{
    return fSi_substrat_thick;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetDetec_shield(G4bool value)
{
    fDetec_shield = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool DetectorConstruction::GetDetec_shield()
{
    return fDetec_shield;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetDetc_shield_thick(G4double value)
{
    fDetc_shield_thick = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double DetectorConstruction::GetDetc_shield_thick()
{
    return fDetc_shield_thick;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetDetec_SiO2_layers(G4bool value)
{
    fDetec_SiO2_layers = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool DetectorConstruction::GetDetec_SiO2_layers()
{
    return fDetec_SiO2_layers;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetSiO2_layers_think(G4double value)
{
    fSiO2_layers_think = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double DetectorConstruction::GetSiO2_layers_think()
{
    return fSiO2_layers_think;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::Setz_lens(G4double value)
{
    fz_lens = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double DetectorConstruction::Getz_lens()
{
    return fz_lens;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::Setz_colr_fltr(G4double value)
{
    colr_fltr = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double DetectorConstruction::Getz_colr_fltr()
{
    return colr_fltr;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::Set_thick_shiel_down(G4double value)
{
    fthick_shiel_down = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double DetectorConstruction::Get_thick_shiel_down()
{
    return fthick_shiel_down;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::Set_thick_shiel_up(G4double value)
{
    fthick_shiel_up = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double DetectorConstruction::Get_thick_shiel_up()
{
    return fthick_shiel_up;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::ConstructSDandField()
{
    G4String filterName, particleName;

    G4SDParticleFilter* alphaFilter =
        new G4SDParticleFilter(filterName="alphaFilter",particleName="alpha");
    G4SDParticleFilter* gammaFilter =
        new G4SDParticleFilter(filterName="gammaFilter",particleName="gamma");
    G4SDParticleFilter* electronFilter =
        new G4SDParticleFilter(filterName="electronFilter",particleName="e-");
    G4SDParticleFilter* positronFilter =
        new G4SDParticleFilter(filterName="positronFilter",particleName="e+");
    G4SDParticleFilter* epFilter = new G4SDParticleFilter(filterName="epFilter");
    epFilter->add(particleName="e-");
    epFilter->add(particleName="e+");

    G4SDParticleFilter* muonminusFilter =
        new G4SDParticleFilter(filterName="muonminusFilter",particleName="mu-");
    G4SDParticleFilter* muonplusFilter =
        new G4SDParticleFilter(filterName="muonplusFilter",particleName="mu+");
    G4SDParticleFilter* muonFilter = new G4SDParticleFilter(filterName="muonFilter");
    muonFilter->add(particleName="mu-");
    muonFilter->add(particleName="mu+");

    G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

    // declare crystal as a MultiFunctionalDetector scorer
    //

    G4MultiFunctionalDetector* cryst = new G4MultiFunctionalDetector("crystal");
    G4SDManager::GetSDMpointer()->AddNewDetector(cryst);

    G4VSensitiveDetector* energy_SD = new EnergySD("EnergySD");
    G4SDManager::GetSDMpointer()->AddNewDetector(energy_SD);
    //logicCryst->SetSensitiveDetector(energy_SD);


    G4VPrimitiveScorer* primitive = new G4PSEnergyDeposit("edep");
    cryst->RegisterPrimitive(primitive);

    primitive = new G4PSCellCharge("chargedep");
    cryst->RegisterPrimitive(primitive);

    primitive = new G4PSTrackLength("TrackLength");  //sum of step lengths of the particles inside the cell.
    cryst->RegisterPrimitive(primitive);

    primitive = new G4PSPassageTrackLength("PassTrackL"); //tracks which pass through the volume are taken into account.
    cryst->RegisterPrimitive(primitive);

    primitive = new G4PSNofStep("nStep");
    cryst->RegisterPrimitive(primitive);




    primitive = new G4PSMinKinEAtGeneration("minEkin");
    cryst->RegisterPrimitive(primitive);

    primitive = new G4PSNofSecondary("nAlpha");
    primitive->SetFilter(alphaFilter);
    cryst->RegisterPrimitive(primitive);

    primitive = new G4PSNofSecondary("nGamma");
    primitive->SetFilter(gammaFilter);
    cryst->RegisterPrimitive(primitive);

    primitive = new G4PSNofSecondary("nElectron");
    primitive->SetFilter(electronFilter);
    cryst->RegisterPrimitive(primitive);

    primitive = new G4PSNofSecondary("nPositron");
    primitive->SetFilter(positronFilter);
    cryst->RegisterPrimitive(primitive);

    primitive = new G4PSNofSecondary("nMuonMinus");
    primitive->SetFilter(muonminusFilter);
    cryst->RegisterPrimitive(primitive);

    primitive = new G4PSNofSecondary("nMuonPlus");
    primitive->SetFilter(muonplusFilter);
    cryst->RegisterPrimitive(primitive);


    primitive = new G4PSEnergyDeposit("edepAlpha");
    primitive->SetFilter(alphaFilter);
    cryst->RegisterPrimitive(primitive);

    primitive = new G4PSEnergyDeposit("edepGamma");
    primitive->SetFilter(gammaFilter);
    cryst->RegisterPrimitive(primitive);

    primitive = new G4PSEnergyDeposit("edepElectron");
    primitive->SetFilter(electronFilter);
    cryst->RegisterPrimitive(primitive);

    primitive = new G4PSEnergyDeposit("edepPositron");
    primitive->SetFilter(positronFilter);
    cryst->RegisterPrimitive(primitive);

    primitive = new G4PSEnergyDeposit("edepMuonMinus");
    primitive->SetFilter(muonminusFilter);
    cryst->RegisterPrimitive(primitive);

    primitive = new G4PSEnergyDeposit("edepMuonPlus");
    primitive->SetFilter(muonplusFilter);
    cryst->RegisterPrimitive(primitive);

    /*
        primitive = new G4PSMinKinEAtGeneration("minEkinAlpha");
        primitive->SetFilter(alphaFilter);
        cryst->RegisterPrimitive(primitive);

        primitive = new G4PSMinKinEAtGeneration("minEkinGamma");
        primitive->SetFilter(gammaFilter);
        cryst->RegisterPrimitive(primitive);

        primitive = new G4PSMinKinEAtGeneration("minEkinElectron");
        primitive->SetFilter(electronFilter);
        cryst->RegisterPrimitive(primitive);

        primitive = new G4PSMinKinEAtGeneration("minEkinPositron");
        primitive->SetFilter(positronFilter);
        cryst->RegisterPrimitive(primitive);

        primitive = new G4PSMinKinEAtGeneration("minEkinMuonMinus");
        primitive->SetFilter(muonminusFilter);
        cryst->RegisterPrimitive(primitive);

        primitive = new G4PSMinKinEAtGeneration("minEkinMuonPlus");
        primitive->SetFilter(muonplusFilter);
        cryst->RegisterPrimitive(primitive);
    */

    //  primitive = new G4PSNofStep("nStep");
    // primitive->SetFilter(epFilter);
    // cryst->RegisterPrimitive(primitive);


    G4PSEnergyDeposit* energyDeposit = new G4PSEnergyDeposit("EnergyDeposit");
    cryst->RegisterPrimitive(energyDeposit);

    SetSensitiveDetector("CrystalLV",energy_SD);
    SetSensitiveDetector("CrystalLV",cryst);

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
