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
// $Id: XFCTSteppingAction.cc 74483 2013-10-09 13:37:06Z gcosmo $
//
/// \file XFCTSteppingAction.cc
/// \brief Implementation of the XFCTSteppingAction class

#include "XFCTSteppingAction.hh"
#include "XFCTEventAction.hh"
#include "XFCTDetectorConstruction.hh"
#include "XFCTHistoManager.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

XFCTSteppingAction::XFCTSteppingAction(XFCTEventAction* eventAction)
: G4UserSteppingAction(),
  fEventAction(eventAction),
  fScoringVolume(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

XFCTSteppingAction::~XFCTSteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void XFCTSteppingAction::UserSteppingAction(const G4Step* step)
{
  if (!fScoringVolume) { 
    const XFCTDetectorConstruction* detectorConstruction
      = static_cast<const XFCTDetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fScoringVolume = detectorConstruction->GetScoringVolume();   
  }

  // get volume of the current step
  G4LogicalVolume* volume 
    = step->GetPreStepPoint()->GetTouchableHandle()
      ->GetVolume()->GetLogicalVolume();
      
  // check if we are in scoring volume
  if (volume != fScoringVolume) return;

  // collect energy deposited in this step
  G4double edepStep = step->GetTotalEnergyDeposit();
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  G4int id = 1;
  
  //if ((edepStep/keV < 26.7) || (27.0 < edepStep/keV < 31.7) || (edepStep/keV > 32.0)) 
  if (edepStep/keV < 26. || edepStep/keV > 32. ) {
  fEventAction->AddEdep(edepStep);
  //plot final state
  analysisManager->FillH1(id,edepStep);
  }

  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

