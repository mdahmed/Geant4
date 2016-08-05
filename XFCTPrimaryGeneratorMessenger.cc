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
// $Id: XFCTPrimarygeneratorMessenger.cc
// GEANT4 tag $Name: xray_fluo-V03-02-00
//
// Author: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//
// History:
// -----------
// 28 Nov 2001 Elena Guardincerri     Created
//
// -------------------------------------------------------------------

#include "XFCTPrimaryGeneratorMessenger.hh"
#include "XFCTPrimaryGeneratorAction.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XFCTPrimaryGeneratorMessenger::XFCTPrimaryGeneratorMessenger(XFCTPrimaryGeneratorAction* XFCTGun)
  :XFCTAction(XFCTGun)
{ 
  /*spectrum = new G4UIcmdWithAString("/gun/spectrum",this);
  spectrum->SetGuidance("Shoot the incident particle with a certain energy spectrum.");
  spectrum->SetGuidance("  Choice : on(default), off");
  spectrum->SetParameterName("choice",true);
  spectrum->SetDefaultValue("on");
  spectrum->SetCandidates("on off");
  spectrum->AvailableForStates(G4State_PreInit,G4State_Idle);*/

  coneBeamCmd = new G4UIcmdWithAString("/gun/coneBeam",this);
  coneBeamCmd->SetGuidance("Choose point source with zero divergence or Cone beam");
  coneBeamCmd->SetGuidance("  Choice : on(default), off");
  coneBeamCmd->SetParameterName("choice",true);
  coneBeamCmd->SetDefaultValue("on");
  coneBeamCmd->SetCandidates("on off");
  coneBeamCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  /*loadPahseSpace=new G4UIcmdWithAString("/gun/loadGunData",this);
  loadPahseSpace->SetGuidance("Load emission from samples form previous runs");
  loadPahseSpace->SetGuidance("Please enter the filename");
  loadPahseSpace->SetParameterName("choice",true);
  loadPahseSpace->AvailableForStates(G4State_Idle);

  loadRayleighData=new G4UIcmdWithABool("/gun/loadRayleighFlag",this);
  loadRayleighData->SetGuidance("Select if data form rayleigh scattering must be loaded");
  loadRayleighData->SetGuidance("To be used before and togheter with /gun/loadGunData");
  loadRayleighData->SetParameterName("Rayleigh Flag",true);
  loadRayleighData->SetDefaultValue(true);
  loadRayleighData->AvailableForStates(G4State_Idle);*/

  G4cout << "XFCTPrimaryGeneratorMessenger created" << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XFCTPrimaryGeneratorMessenger::~XFCTPrimaryGeneratorMessenger()
{
  //delete spectrum;
  delete coneBeamCmd;
  //delete loadPahseSpace;
  //delete loadRayleighData;
  G4cout << "XFCTPrimaryGeneratorMessenger deleted" << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XFCTPrimaryGeneratorMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{ 
  
 /*if( command == spectrum )
   { XFCTAction->SetSpectrum(newValue);}*/ 
 if( command == coneBeamCmd )
   { XFCTAction->SetConeBeam(newValue);}
 /*if( command == loadPahseSpace )
   { XFCTAction->ActivatePhaseSpace(newValue);}*/
 /*if( command == loadRayleighData )
   { 
     G4cout << "newValue: " << newValue << G4endl;

     G4bool newRayFlag = loadRayleighData->GetNewBoolValue(newValue);
     XFCTAction->SetRayleighFlag(newRayFlag);
   }*/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

