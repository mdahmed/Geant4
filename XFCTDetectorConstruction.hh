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
// $Id: XFCTDetectorConstruction.hh 69565 2013-05-08 12:35:31Z gcosmo $
//
/// \file XFCTDetectorConstruction.hh
/// \brief Definition of the XFCTDetectorConstruction class

#ifndef XFCTDetectorConstruction_h
#define XFCTDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;
class XFCTDetectorMessenger;

/// Detector construction class to define materials and geometry.

class XFCTDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    XFCTDetectorConstruction();
    ~XFCTDetectorConstruction();

    //void DefineMaterials();
    G4VPhysicalVolume* DefineVolumes();
    virtual G4VPhysicalVolume* Construct();
    
    void			SetFilter1Material (G4String );
    void			SetFilter2Material (G4String );
    
    void			SetFilterToSampleDistance(G4double value);
    void			SetDetToSampleDistance(G4double value);
    void			SetFilter1Thickness(G4double value);
    void			SetFilter2Thickness(G4double value);
    
    void			SetNanoParticleConc(G4double value);
    
    void			SetCollimRadius(G4double value);
    void			SetCollimThickness(G4double value);
    
    //void			PrintParameters();

    G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }
    
  private:
    G4LogicalVolume*   		flogicFilter1;     // pointer to the logical Filter1
    G4LogicalVolume*   		flogicFilter2;     // pointer to the logical Filter2

    G4Material*				fMatFilter1;  // pointer to the Filter1  material
    G4Material*				fMatFilter2;  // pointer to the Filter2  material

	G4double				fFilterToSampleDis;
	G4double            	fDetToSampleDis;
	G4double            	fThicknessFilter1;
	G4double            	fThicknessFilter2;
	
	G4double            	fGoldFmass;
	
	G4double            	fCollimRadius;
	G4double            	fCollimThickness;
	
	//G4ThreeVector			fdetPosition;
	//G4Transform3D			fdetTransform;
    XFCTDetectorMessenger*  fMessenger;   // messenger


    G4bool  fCheckOverlaps;

  protected:
    G4LogicalVolume*  fScoringVolume;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

