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
// $Id: XFCTDetectorConstruction.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file XFCTDetectorConstruction.cc
/// \brief Implementation of the XFCTDetectroConstruction class

#include "XFCTDetectorConstruction.hh"
#include "XFCTDetectorMessenger.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4SystemOfUnits.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

//--------------------------------------------------------------------
//								Varibles
//--------------------------------------------------------------------
/*	Components		Quantity		Unit Catagory		Status
 *	----------------------------------------------------------------		
 * 	Filter 1		Material		Material			Variable
 * 					Thickness		Length				Variable
 * 	----------------------------------------------------------------
 *	Filter 2		Material		Material			Variable
 * 					Thickness		Length				Variable
 * 	----------------------------------------------------------------
 * 	Au Conc.		Frac. mass		Value				Variable(?)
 * 	----------------------------------------------------------------
 * 	Collimator		Radius			Length				Variable
 * 					Thickness		Length				Variable
 * 	----------------------------------------------------------------
 * 	Detector		Material		Sensitivity			Fixed
 * 					Distance		Length				Variable
 * 	----------------------------------------------------------------
 * 
 */


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


XFCTDetectorConstruction::XFCTDetectorConstruction()
: G4VUserDetectorConstruction(),
	//fdetPosition(0),
	//fdetTransform(0),
	fScoringVolume(0)
{
  fFilterToSampleDis = 15.0*cm;
  fDetToSampleDis = 10.*cm;
  fThicknessFilter1 = 0.2*cm;
  fThicknessFilter2 = 0.2*cm;
  fGoldFmass = 0.01;
  fCollimRadius = 0.25*cm;
  fCollimThickness = 4.*cm;
  fMessenger = new XFCTDetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

XFCTDetectorConstruction::~XFCTDetectorConstruction()
{
	delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* XFCTDetectorConstruction::Construct()
{
  // Define volumes
  return DefineVolumes();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume* XFCTDetectorConstruction::DefineVolumes()
{
//---------------------------------------------------------
// 					Cleanup old geometry
//---------------------------------------------------------  
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
  
//----------------------------------------------------------------------
//							Define Materials
//----------------------------------------------------------------------

 G4NistManager* nist = G4NistManager::Instance();
  
  G4Material* LabAir = nist->FindOrBuildMaterial("G4_AIR");
  
  G4double density, fractionmass, ncomponents;
  
  G4Material* H2O  = nist->FindOrBuildMaterial("G4_WATER");
  G4Material* Gold  = nist->FindOrBuildMaterial("G4_Au");
  
  G4double fmass_H2O = 1. - fGoldFmass;
  
  // Direct density input from NIST
  G4double nistDensity_Gold = 19.32*g/cm3;
  G4double nistDensity_H2O = 1.*g/cm3;
  
  density = fGoldFmass*nistDensity_Gold + fmass_H2O*nistDensity_H2O;
  G4Material* GoldSol = new G4Material("GoldSolution", density,2);
  GoldSol->AddMaterial(H2O, fractionmass=fmass_H2O);
  GoldSol->AddMaterial(Gold, fractionmass=fGoldFmass);

  
  G4bool isotopes = false;
    G4Element* Fe = nist->FindOrBuildElement("Fe" , isotopes); 
	G4Element* Ni = nist->FindOrBuildElement("Ni", isotopes);
	G4Element* Co = nist->FindOrBuildElement("Co", isotopes);
	G4Element* C  =  nist->FindOrBuildElement("C" , isotopes); 
	G4Element* Si = nist->FindOrBuildElement("Si", isotopes);
	G4Element* Mn = nist->FindOrBuildElement("Mn", isotopes);
	G4Element* Cr = nist->FindOrBuildElement("Cr" , isotopes);
	
  // Detector head made of Kovar steel
  
	G4Material* KobarSteel = new G4Material("KovarSteel",
	 density= 8.0*g/cm3, ncomponents=6);
	KobarSteel->AddElement(Fe, fractionmass=0.5349);
	KobarSteel->AddElement(Ni, fractionmass=0.2900);
	KobarSteel->AddElement(Co, fractionmass=0.1700);
	KobarSteel->AddElement(C,  fractionmass=0.0001);
	KobarSteel->AddElement(Si, fractionmass=0.0020);
	KobarSteel->AddElement(Mn, fractionmass=0.0030);
	
	
	
	G4Material* StainlessSteel = new G4Material("StainlessSteel",
	 density= 8.06*g/cm3, ncomponents=6);
	StainlessSteel->AddElement(C,  fractionmass=0.001);
	StainlessSteel->AddElement(Si, fractionmass=0.007);
	StainlessSteel->AddElement(Cr, fractionmass=0.180);
	StainlessSteel->AddElement(Mn, fractionmass=0.010);
	StainlessSteel->AddElement(Fe, fractionmass=0.712);
	StainlessSteel->AddElement(Ni, fractionmass=0.090);
	
	
	fMatFilter1 = nist->FindOrBuildMaterial("G4_AIR");
	fMatFilter2 = nist->FindOrBuildMaterial("G4_AIR");
	
	G4Material* Lead = nist->FindOrBuildMaterial("G4_Pb");
	
	G4Material* crys_mat = nist->FindOrBuildMaterial("G4_CADMIUM_TELLURIDE");
	
	// All known materials
	G4Material* detShield_mat = nist->FindOrBuildMaterial("G4_Pb");
	G4Material* detHead_mat = nist->FindOrBuildMaterial("G4_Ni");
	G4Material* detWindow_mat = nist->FindOrBuildMaterial("G4_Be");
	
	G4Material* sampleHolder_mat = nist->FindOrBuildMaterial("G4_POLYETHYLENE");
	
   
// Option to switch on/off checking of volumes overlaps
//
  G4bool checkOverlaps = true;

//---------------------------------------------------------
// 							World
//---------------------------------------------------------  
  G4double world_sizeXY = 200.*cm;
  G4double world_sizeZ  = 200.*cm;
  
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        LabAir,           //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking
                     
 
//---------------------------------------------------------
// 					Cone beam generator
//---------------------------------------------------------    
        
  // Cone filter aperture       
  G4double conFilAir_rmina =  0.*cm, conFilAir_rmaxa = 0.9*0.5*cm;
  G4double conFilAir_rminb =  0.*cm, conFilAir_rmaxb = 2.*0.5*cm;
  G4double conFilter_hz = 5.1*cm;
  G4double conFilAir_phimin = 0.*deg, conFilAir_phimax = 360.*deg;
  G4Cons* conFilAir_solid =    
    new G4Cons("ConeAperture", 
    conFilAir_rmina, conFilAir_rmaxa, conFilAir_rminb, conFilAir_rmaxb, conFilter_hz*0.5,
    conFilAir_phimin, conFilAir_phimax);
                      
  G4LogicalVolume* conFilAir_logic =                         
    new G4LogicalVolume(conFilAir_solid,         //its solid
                        LabAir,          //its material
                        "ConeAperture");           //its name
                        
  G4VisAttributes* conFilAir_color = new G4VisAttributes(G4Colour(0.5,0.5,0.5));
  conFilAir_logic->SetVisAttributes(conFilAir_color);
  //conFilAir_logic->SetVisAttributes (G4VisAttributes::Invisible);
  
  // Cone filter cube      
  G4double conFil_X =  7.*cm, conFil_Y = 7.*cm;  
  
  G4Box* conFil_solid =    
    new G4Box("ConeFilterCube",                    //its name
        0.5*conFil_X, 0.5*conFil_Y, 0.5*conFilter_hz); //its size
    
                      
  G4LogicalVolume* conFil_logic =                         
    new G4LogicalVolume(conFil_solid,         //its solid
                        Lead,          //its material
                        "ConeFilterCube");           //its name
                        
  G4VisAttributes* conFil_color = new G4VisAttributes(G4Colour(0.5,0.5,0.5));
  conFil_logic->SetVisAttributes(conFil_color);
  
  // Place cone filter aperture inside cone filter cube
               
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0, 0.,0.), //at position
                    conFilAir_logic,         //its logical volume
                    "Filter-aperture",            //its name
                    conFil_logic,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
                    
  // Place cone filter inside world
  
  G4double conFil_zDistance = 0.*cm;
  G4ThreeVector conFil_position = 
  G4ThreeVector(0, 0., conFil_zDistance + 0.5*conFilter_hz);
                    
  new G4PVPlacement(0,                       //no rotation
                    conFil_position, //at position
                    conFil_logic,         //its logical volume
                    "ConeFilter",            //its name
                    logicWorld,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
                    
          
//---------------------------------------------------------
// 						Phase Space
//---------------------------------------------------------  
  
       
  G4double phaseSpace_hz = 1.0e-4*cm; 
 
  G4Tubs* phaseSpace_solid = 
	 new G4Tubs("PhaseSpace_solid", 0., 1.5*fDetToSampleDis,
	  phaseSpace_hz, 0., 360.*deg);
	 

	
  G4LogicalVolume* phaseSpace_logic =                         
    new G4LogicalVolume(phaseSpace_solid,      //its solid
                        LabAir,        //its material
                        "PhaseSpace_logic");   //its name
                        
  G4VisAttributes* phaseSpace_color = new G4VisAttributes(G4Colour(0.,0.,1.));
  phaseSpace_color->SetVisibility(false);
  phaseSpace_logic->SetVisAttributes(phaseSpace_color);
  
  G4double phaseSpace_zDistance = conFil_zDistance + 
  conFilter_hz + 1.0e-4*cm + phaseSpace_hz*0.5;
  G4ThreeVector phaseSpace_pos = G4ThreeVector(0, 0., phaseSpace_zDistance);
                    
  new G4PVPlacement(0,                      //no rotation
                    phaseSpace_pos, 			//at position
                    phaseSpace_logic,         	//its logical volume
                    "PhaseSpace",            	//its name
                    logicWorld,               //its mother  volume
                    false,                  //no boolean operation
                    0,                      //copy number
                    checkOverlaps);         //overlaps checking
                    
                    
//---------------------------------------------------------
// 							Filter
//---------------------------------------------------------  
  
  // Make x-y dimensions of all the filters same as cone filter
  
  //....oooOO0OOooo........oooOO0OOooo.......
  // Filter 1
  
  //G4double fThicknessFilter1 = 0.2*cm; // should be variable      
  
  G4Box* filter1_solid =    
    new G4Box("Filter1_solid",                    //its name
        0.5*conFil_X, 0.5*conFil_Y, 0.5*fThicknessFilter1); //its size
	
  flogicFilter1 =                         
    new G4LogicalVolume(filter1_solid,      //its solid
                        fMatFilter1,        //its material
                        "Filter1_logic");   //its name
                        
  G4VisAttributes* filter1_color = new G4VisAttributes(G4Colour(1.,0.,0.));
  flogicFilter1->SetVisAttributes(filter1_color);
  
  G4double filter1_zDistance = phaseSpace_zDistance + 0.5*phaseSpace_hz
    + 1.0e-4*cm + 0.5*fThicknessFilter1;
  G4ThreeVector filter1_pos = G4ThreeVector(0, 0., filter1_zDistance);
                    
  new G4PVPlacement(0,                      //no rotation
                    filter1_pos, 			//at position
                    flogicFilter1,         	//its logical volume
                    "Filter_1",            	//its name
                    logicWorld,               //its mother  volume
                    false,                  //no boolean operation
                    0,                      //copy number
                    checkOverlaps);         //overlaps checking
                    
  //....oooOO0OOooo........oooOO0OOooo.......
  // Filter 2
  
  G4Box* filter2_solid =    
    new G4Box("Filter2_solid",                    //its name
        0.5*conFil_X, 0.5*conFil_Y, 0.5*fThicknessFilter2); //its size
	
  flogicFilter2 =                         
    new G4LogicalVolume(filter2_solid,      //its solid
                        fMatFilter2,        //its material
                        "Filter2_logic");         //its name
                        
  G4VisAttributes* filter2_color = new G4VisAttributes(G4Colour(0.,1.,0.));
  flogicFilter2->SetVisAttributes(filter2_color);
  
  G4double filter2_zDistance = filter1_zDistance + 0.5*fThicknessFilter1
    + 1.0e-4*cm + 0.5*fThicknessFilter2;

  G4ThreeVector filter2_pos = G4ThreeVector(0, 0., filter2_zDistance);
                    
  new G4PVPlacement(0,                      //no rotation
                    filter2_pos, 			//at position
                    flogicFilter2,         	//its logical volume
                    "Filter_2",            	//its name
                    logicWorld,               //its mother  volume
                    false,                  //no boolean operation
                    0,                      //copy number
                    checkOverlaps);         //overlaps checking
  
                    
                    
//---------------------------------------------------------
// 				Sample (nano-particle solution)
//---------------------------------------------------------    
	
	// Liquid nano-particle sample       
	G4double nanoSol_Ri =  0.*cm;
	G4double nanoSol_Ro =  0.5*cm;
	G4double nanoSol_hz = 2.0*cm;
	G4double nanoSol_startAngle = 0.*deg;
	G4double nanoSol_spanningAngle = 360.*deg;
	
	G4Tubs* solidnanoSol = 
	 new G4Tubs("Sample_Soln", nanoSol_Ri, nanoSol_Ro,
	 nanoSol_hz, nanoSol_startAngle, nanoSol_spanningAngle);
	
	G4LogicalVolume* logicnanoSol = 
	 new G4LogicalVolume(solidnanoSol,         	//its solid
	 GoldSol,          						//its material
	"Sample_Soln");        						//its name
	
	G4VisAttributes* solnColor = new G4VisAttributes(G4Colour(1.,0.,0.));
	 logicnanoSol->SetVisAttributes(solnColor);
	
	// Sample holder! build with refernce to sample!
	
	G4double SmplHolder_Ro = nanoSol_Ro + 0.1*cm;
	G4Tubs* solidholder = 
	 new G4Tubs("Sample_Holder", nanoSol_Ri, SmplHolder_Ro,
	 nanoSol_hz, nanoSol_startAngle, nanoSol_spanningAngle);
	
	G4LogicalVolume* logicholder = 
	 new G4LogicalVolume(solidholder,         	//its solid
	 sampleHolder_mat,          						//its material
	"Sample_Holder");        						//its name
	
	// Place sample inside holder
	
	new G4PVPlacement(0,						//no rotation
					G4ThreeVector(0, 0, 0.*cm), //at position
					logicnanoSol,             	//its logical volume
					"Sample_Tube",                	//its name
					logicholder,                	//its mother volume
					false,                   	//no boolean operation
					0,                       	//copy number
					checkOverlaps);          	//overlaps checking
	
	G4double sample_phi = 0.*deg;
    G4RotationMatrix sample_rotm  = G4RotationMatrix();
    sample_rotm.rotateX(90.*deg); 
    sample_rotm.rotateZ(sample_phi);
        
    G4double SourceToSampleDis = filter2_zDistance + 0.5*fThicknessFilter2
    + 0.5*nanoSol_hz + fFilterToSampleDis;
    
    G4ThreeVector sample_position = G4ThreeVector(0., 0., SourceToSampleDis);
    G4Transform3D sample_transform = G4Transform3D(sample_rotm,sample_position);
	
	new G4PVPlacement(sample_transform,		// its rotation and position
					logicholder,           //its logical volume
					"Sample",                //its name
					logicWorld,                //its mother  volume
					false,                   //no boolean operation
					0,                       //copy number
					checkOverlaps);          //overlaps checking
					
					
					
//---------------------------------------------------------
// 							Detector
//---------------------------------------------------------  
	
	// mu -> thickness, L -> length, w -> width, R -> Radius
	// Ri -> innerRadius, Ro -> outerRadius
	// SPhi -> startAngle, DPhi -> spanningAngle
	
	G4double detShield_mu = 0.2*cm;  		// From Nivedh's diss
	G4double detBase_mu = 0.15*cm;	// From Nivedh's diss
	G4double detHead_mu = 0.025*cm;  // From Nivedh's diss
	
	G4double detBase_L = 3.81*cm; 	// From Amptek
	G4double detHead_L = 0.86*cm; 	// From Amptek
	G4double detWindow_L = 0.01*cm; // From Amptek
	
	G4double crys_w = 0.3*cm; // From Amptek
	G4double crys_L = 0.1*cm; // From Amptek
	
	G4double wind_to_crys_L 
						= (0.127 - detWindow_L*0.5)*cm; // From Amptek
	
	G4double detBase_R = 0.89*cm; // From Amptek
	G4double detHead_R = 0.70*cm; // From Amptek
	
	//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........
	// Detector shield (Pb shielding): Mother of entire detector
	
	G4double detShield_Ri = 0.*cm; 	// Solid
	G4double detShield_Ro = detBase_R + detShield_mu + 0.1*cm;
	G4double detShield_h = detBase_L + detHead_L + 0.1*cm; // Extra room
	G4double detShield_SPhi = 0.*deg;
	G4double detShield_DPhi = 360.*deg;
	
	G4Tubs* detShield_solid = new G4Tubs("Shield_solid", 
							detShield_Ri, detShield_Ro, detShield_h*0.5, 
							detShield_SPhi, detShield_DPhi);
	
	G4LogicalVolume* shield_logic = new G4LogicalVolume(
			detShield_solid, detShield_mat, "Shield_logic");
	
	G4VisAttributes* shield_color = 
							new G4VisAttributes(G4Colour(0.5,0.5,0.5));
	 shield_logic->SetVisAttributes(shield_color);
	 
	
	// Chamber air inside shield
	
	G4double chamAir_Ri = 0.*cm; 	// Solid
	G4double chamAir_Ro = detBase_R + 0.1*cm;
	G4double chamAir_h = detShield_h;
	G4double chamAir_SPhi = 0.*deg;
	G4double chamAir_DPhi = 360.*deg;
	
	G4Tubs* chamAir_solid = new G4Tubs("Chamber_solid", 
							chamAir_Ri, chamAir_Ro, chamAir_h*0.5, 
							chamAir_SPhi, chamAir_DPhi);
	
	G4LogicalVolume* chamAir_logic = new G4LogicalVolume(
			chamAir_solid, LabAir, "Chamber_logic");
	
	G4VisAttributes* chamAir_color = new G4VisAttributes(G4Colour(0.,0.,0.));
	chamAir_logic->SetVisAttributes(chamAir_color);
	 
	//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........
	// Detector base
	
	G4double detBase_Ri 	= 0.*cm; 	// Solid
	G4double detBase_Ro 	= detBase_R;
	G4double detBase_h 		= detBase_L; // Extra room
	G4double detBase_SPhi 	= 0.*deg;
	G4double detBase_DPhi 	= 360.*deg;
	
	G4Tubs* detBase_solid = new G4Tubs("Base_solid", 
							detBase_Ri, detBase_Ro, detBase_h*0.5, 
							detBase_SPhi, detBase_DPhi);
	
	G4LogicalVolume* base_logic = new G4LogicalVolume(
			detBase_solid, KobarSteel, "Base_logic");
	
	G4VisAttributes* base_color = 
							new G4VisAttributes(G4Colour(0.,0.,1.0));
	 base_logic->SetVisAttributes(base_color);
	 
	 // Air inside detector base
	
	G4double baseAir_Ri 	= 0.*cm; 	// Solid
	G4double baseAir_Ro 	= detBase_R - detBase_mu;
	G4double baseAir_h 		= detBase_L; // Extra room
	G4double baseAir_SPhi 	= 0.*deg;
	G4double baseAir_DPhi 	= 360.*deg;
	
	G4Tubs* baseAir_solid = new G4Tubs("BaseAir_solid", 
							baseAir_Ri, baseAir_Ro, baseAir_h*0.5, 
							baseAir_SPhi, baseAir_DPhi);
	
	G4LogicalVolume* baseAir_logic = new G4LogicalVolume(
			baseAir_solid, LabAir, "BaseAir_logic");
	
	G4VisAttributes* baseAir_color = 
							new G4VisAttributes(G4Colour(0.,0.,0.));
	 baseAir_logic->SetVisAttributes(baseAir_color);
	 
	 //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........
	 // Detector head
	
	G4double detHead_Ri 	= 0.*cm; 	// Solid
	G4double detHead_Ro 	= detHead_R;
	G4double detHead_h 		= detHead_L; // Extra room
	G4double detHead_SPhi 	= 0.*deg;
	G4double detHead_DPhi 	= 360.*deg;
	
	G4Tubs* detHead_solid = new G4Tubs("Head_solid", 
							detHead_Ri, detHead_Ro, detHead_h*0.5, 
							detHead_SPhi, detHead_DPhi);
	
	G4LogicalVolume* head_logic = new G4LogicalVolume(
			detHead_solid, detHead_mat, "Head_logic");
	
	G4VisAttributes* head_color = 
							new G4VisAttributes(G4Colour(0.,1.0,0.));
	 head_logic->SetVisAttributes(head_color);
	 
	// Air inside detector head
	
	G4double headAir_Ri 	= 0.*cm; 	// Solid
	G4double headAir_Ro 	= detHead_R - detHead_mu;
	G4double headAir_h 		= detHead_L; // Extra room
	G4double headAir_SPhi 	= 0.*deg;
	G4double headAir_DPhi 	= 360.*deg;
	
	G4Tubs* headAir_solid = new G4Tubs("HeadAir_solid", 
							headAir_Ri, headAir_Ro, headAir_h*0.5, 
							headAir_SPhi, headAir_DPhi);
	
	G4LogicalVolume* headAir_logic = new G4LogicalVolume(
			headAir_solid, LabAir, "HeadAir_logic");
	
	G4VisAttributes* headAir_color = 
							new G4VisAttributes(G4Colour(0.,0.,0.));
	 headAir_logic->SetVisAttributes(headAir_color);
	 
	 // Crystal
	
	G4Box* crystal_solid = 
	 new G4Box("Crystal", 0.5*crys_w, 0.5*crys_w, 0.5*crys_L);
	
	G4LogicalVolume* crystal_logic = 
				new G4LogicalVolume(crystal_solid, crys_mat, "Crystal");           						//its name
	
	G4VisAttributes* crys_color = 
							new G4VisAttributes(G4Colour(1.,0,0));
	 crystal_logic->SetVisAttributes(crys_color);
	 
	//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........
	// Detector window (Be)
	
	G4double detWindow_Ri = 0; 			//Solid
	G4double detWindow_Ro = detHead_Ro;	//fixed
	G4double detWindow_h = detWindow_L;				//fixed values
	G4double detWindow_SPhi = 0.*deg;			//fixed values
	G4double detWindow_DPhi = 360.*deg;	//fixed values
	
	G4Tubs* detWindow_solid = new G4Tubs("DetWindow_solid", 
				detWindow_Ri,
				detWindow_Ro,
				detWindow_h*0.5,
				detWindow_SPhi,
				detWindow_DPhi
				);
	
	G4LogicalVolume* detWindow_logic = new G4LogicalVolume(
			detWindow_solid,         				//its solid
			detWindow_mat,          				//its material
			"DetWindow_logic");        			//its name
	
	G4VisAttributes* detWindow_color = new G4VisAttributes(G4Colour(1.,0.,0.));
	 detWindow_logic->SetVisAttributes(detWindow_color);
	 
	 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......	
	 
	//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........
	// window: first daughter of chamber air mother volume
	
	G4ThreeVector detwindow_pos = 
				G4ThreeVector(0, 0, -(detShield_h - detWindow_h)*0.5);
	new G4PVPlacement(0,						//0 deg rotation
					detwindow_pos,        			//at position
					detWindow_logic,            		//its logical volume
					"window-chamber",            //its name
					chamAir_logic,               //its mother  volume
					false,                   	//no boolean operation
					0,                       	//copy number
					checkOverlaps);          	//overlaps checking
	
	//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........
	// Head-with-crystal: second daughter of chamber air mother volume
	
	// Place crystal in the air inside head
	
	G4ThreeVector crystal_pos = 
				G4ThreeVector(0, 0, -detHead_L*0.5 + wind_to_crys_L);
	new G4PVPlacement(0,						//0 deg rotation
					crystal_pos,        			//at position
					crystal_logic,            		//its logical volume
					"Crystal-headAir",            //its name
					headAir_logic,               //its mother  volume
					false,                   	//no boolean operation
					0,                       	//copy number
					checkOverlaps);          	//overlaps checking
	
	// Place crystal-air inside head
	
	G4ThreeVector headAir_pos = G4ThreeVector(0, 0, 0);
	new G4PVPlacement(0,						//0 deg rotation
					headAir_pos,        			//at position
					headAir_logic,            		//its logical volume
					"Crystal-Air-Head",            //its name
					head_logic,               //its mother  volume
					false,                   	//no boolean operation
					0,                       	//copy number
					checkOverlaps);          	//overlaps checking
					
	// Place head with crystal inside chamber air
	G4ThreeVector head_pos = 
				G4ThreeVector(0, 0, -(detShield_h - detHead_h)*0.5 + 0.1*cm);
	new G4PVPlacement(0,						//0 deg rotation
					head_pos,        			//at position
					head_logic,            		//its logical volume
					"Head-Crystal-Chamber",            //its name
					chamAir_logic,               //its mother  volume
					false,                   	//no boolean operation
					0,                       	//copy number
					checkOverlaps);          	//overlaps checking
					
	//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........
	// Base: third daughter of chamber air mother volume
	
	// Place air inside base
	G4ThreeVector baseAir_pos = G4ThreeVector(0, 0, 0);
	new G4PVPlacement(0,						//0 deg rotation
					baseAir_pos,        			//at position
					baseAir_logic,            		//its logical volume
					"Base-Air",            //its name
					base_logic,               //its mother  volume
					false,                   	//no boolean operation
					0,                       	//copy number
					checkOverlaps);          	//overlaps checking
	
	// Place air-filled base inside chamber air
	G4ThreeVector base_pos = 
				G4ThreeVector(0, 0, (detShield_h - detBase_h)*0.5);
	new G4PVPlacement(0,						//0 deg rotation
					base_pos,        			//at position
					base_logic,            		//its logical volume
					"Base-Shield",            //its name
					chamAir_logic,               //its mother  volume
					false,                   	//no boolean operation
					0,                       	//copy number
					checkOverlaps);          	//overlaps checking
	
	
	// Place chamber air inside shield
	G4ThreeVector chamber_pos = G4ThreeVector(0, 0, 0);
	new G4PVPlacement(0,						//0 deg rotation
					chamber_pos,        		//at position
					chamAir_logic,            	//its logical volume
					"Shield-chamber",           //its name
					shield_logic,               //its mother  volume
					false,                   	//no boolean operation
					0,                       	//copy number
					checkOverlaps);          	//overlaps checking
	
	// Place Shield inside world
	G4ThreeVector fdetPosition = 
				G4ThreeVector(fDetToSampleDis + detShield_h*0.5, 0, SourceToSampleDis);
	
	G4double detector_phi = 0.*deg;
    G4RotationMatrix detector_rotm  = G4RotationMatrix();
    detector_rotm.rotateY(90.*deg); 
    detector_rotm.rotateZ(detector_phi);
    G4Transform3D fdetTransform = G4Transform3D(detector_rotm,fdetPosition);
	
	new G4PVPlacement(fdetTransform,		//90 deg rotation
					//fdetPosition,        		//at position
					shield_logic,            	//its logical volume
					"Detector",              	//its name
					logicWorld,                	//its mother  volume
					false,                   	//no boolean operation
					0,                       	//copy number
					checkOverlaps);          	//overlaps checking
					
					
					
//---------------------------------------------------------
// 						Collimator
//---------------------------------------------------------    
  
  G4double det_to_collimator = 0.2*cm; // should be variable
  //G4double collimator_radius = 0.25*cm; // should be variable
  //G4double collimator_length = 4.*cm;  // should be variable
				 
  G4double collim_XY = 3.*cm;
  
  
	G4double collimApert_Ri = 0.*cm; 	// Solid
	//G4double collimApert_h = collimator_length;
	G4double collimApert_SPhi = 0.*deg;
	G4double collimApert_DPhi = 360.*deg;
	
	
	
	G4Tubs* collimApert_solid = new G4Tubs("CollimApert_solid", 
				collimApert_Ri,
				fCollimRadius,
				fCollimThickness*0.5,
				collimApert_SPhi,
				collimApert_DPhi
				);
	
	G4LogicalVolume* collimApert_logic = new G4LogicalVolume(
			collimApert_solid,         				//its solid
			LabAir,          				//its material
			"CollimApert_logic");       			//its name
	
	G4VisAttributes* collimApert_color = new G4VisAttributes(G4Colour(1.,0.,0.));
	 collimApert_logic->SetVisAttributes(collimApert_color);
  
  
  G4Box* collimBox_solid =    
    new G4Box("CollimatorBox_solid",                    //its name
        0.5*collim_XY, 0.5*collim_XY, 0.5*fCollimThickness); //its size
      
  G4LogicalVolume* collimBox_logic =                         
    new G4LogicalVolume(collimBox_solid,            //its solid
                        StainlessSteel,             //its material
                        "CollimatorBox_logic");         //its name
                        
   G4VisAttributes* collimBox_color = new G4VisAttributes(G4Colour(1.,0.,0.));
	 collimBox_logic->SetVisAttributes(collimBox_color);
                        
  
   new G4PVPlacement(0,							//0 deg rotation
					G4ThreeVector(0., 0., 0.),   //at position
					collimApert_logic,            	//its logical volume
					"CollimBox-CollimApert",    //its name
					collimBox_logic,                	//its mother  volume
					false,                   	//no boolean operation
					0,                       	//copy number
					checkOverlaps);          	//overlaps checking
					
   G4ThreeVector collimator_pos = 
				G4ThreeVector(fDetToSampleDis -
				 det_to_collimator - fCollimThickness*0.5, 0., SourceToSampleDis);
				 
    G4Transform3D collimator_transform = 
							G4Transform3D(detector_rotm,collimator_pos);
				 
   
   new G4PVPlacement(collimator_transform,		//90 deg rotation
					//collimator_pos,        		//at position
					collimBox_logic,            	//its logical volume
					"Collimator",              	//its name
					logicWorld,                	//its mother  volume
					false,                   	//no boolean operation
					0,                       	//copy number
					checkOverlaps);          	//overlaps checking             
                
  // Set Detector as scoring volume
  //
  fScoringVolume = crystal_logic;

  //
  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void XFCTDetectorConstruction::SetFilter1Material(G4String materialName)
{
  G4NistManager* nistManager = G4NistManager::Instance();

  G4Material* pttoMaterial = 
              nistManager->FindOrBuildMaterial(materialName);

  if (fMatFilter1 != pttoMaterial) {
     if ( pttoMaterial ) {
        fMatFilter1 = pttoMaterial;
        if (flogicFilter1) flogicFilter1->SetMaterial(fMatFilter1);
        G4cout 
          << G4endl 
          << "----> Filter 1 is made of " << materialName << G4endl;
     } else {
        G4cout 
          << G4endl 
          << "-->  WARNING from SetFilter1Material : "
          << materialName << " not found" << G4endl;
     }
  }
}

void XFCTDetectorConstruction::SetFilter2Material(G4String materialName)
{
  G4NistManager* nistManager = G4NistManager::Instance();

  G4Material* pttoMaterial = 
              nistManager->FindOrBuildMaterial(materialName);

  if (fMatFilter2 != pttoMaterial) {
     if ( pttoMaterial ) {
        fMatFilter2 = pttoMaterial;
        if (flogicFilter2) flogicFilter2->SetMaterial(fMatFilter2);
        G4cout 
          << G4endl 
          << "----> Filter 2 is made of " << materialName << G4endl;
     } else {
        G4cout 
          << G4endl 
          << "-->  WARNING from SetFilter2Material : "
          << materialName << " not found" << G4endl;
     }
  }
}


void XFCTDetectorConstruction::SetFilter1Thickness(G4double value)
{
  fThicknessFilter1 = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void XFCTDetectorConstruction::SetFilter2Thickness(G4double value)
{
  fThicknessFilter2 = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void XFCTDetectorConstruction::SetFilterToSampleDistance(G4double value)
{
  fFilterToSampleDis = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void XFCTDetectorConstruction::SetDetToSampleDistance(G4double value)
{
  fDetToSampleDis = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}


void XFCTDetectorConstruction::SetNanoParticleConc(G4double value)
{
  fGoldFmass = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void XFCTDetectorConstruction::SetCollimRadius(G4double value)
{
  fCollimRadius = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void XFCTDetectorConstruction::SetCollimThickness(G4double value)
{
  fCollimThickness = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
