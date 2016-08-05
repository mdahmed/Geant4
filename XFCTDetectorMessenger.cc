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
// $Id: B2aDetectorMessenger.cc 69706 2013-05-13 09:12:40Z gcosmo $
// 
/// \file B2aDetectorMessenger.cc
/// \brief Implementation of the B2aDetectorMessenger class

#include "XFCTDetectorMessenger.hh"
#include "XFCTDetectorConstruction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

XFCTDetectorMessenger::XFCTDetectorMessenger(XFCTDetectorConstruction* Det)
 : G4UImessenger(),
   fDetectorConstruction(Det)
{
  fXFCTDirectory = new G4UIdirectory("/XFCT/");
  fXFCTDirectory->SetGuidance("UI commands specific to this example.");

  fDetDirectory = new G4UIdirectory("/XFCT/det/");
  fDetDirectory->SetGuidance("Detector construction control");

  fFilter1MatCmd = new G4UIcmdWithAString("/XFCT/det/setFilter1Material",this);
  fFilter1MatCmd->SetGuidance("Select Material of the Filter 1.");
  fFilter1MatCmd->SetParameterName("choice",false);
  fFilter1MatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fFilter2MatCmd = new G4UIcmdWithAString("/XFCT/det/setFilter2Material",this);
  fFilter2MatCmd->SetGuidance("Select Material of the Filter 2.");
  fFilter2MatCmd->SetParameterName("choice",false);
  fFilter2MatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fThickFil1Cmd = new G4UIcmdWithADoubleAndUnit("/XFCT/det/setFilter1Thickness",this);
  fThickFil1Cmd->SetGuidance("Set thickness of Filter 1");
  fThickFil1Cmd->SetParameterName("Filter1Thickness",false);
  fThickFil1Cmd->SetUnitCategory("Length");
  fThickFil1Cmd->AvailableForStates(G4State_Idle);
  
  fThickFil2Cmd = new G4UIcmdWithADoubleAndUnit("/XFCT/det/setFilter2Thickness",this);
  fThickFil2Cmd->SetGuidance("Set thickness of Filter 2");
  fThickFil2Cmd->SetParameterName("Filter2Thickness",false);
  fThickFil2Cmd->SetUnitCategory("Length");
  fThickFil2Cmd->AvailableForStates(G4State_Idle);
  
  fFilterToSmplCmd = new G4UIcmdWithADoubleAndUnit("/XFCT/det/FilterToSmplDis",this);
  fFilterToSmplCmd->SetGuidance("Set Filter to sample distance");
  fFilterToSmplCmd->SetParameterName("FilToSmplDis",false);
  fFilterToSmplCmd->SetUnitCategory("Length");
  fFilterToSmplCmd->AvailableForStates(G4State_Idle);
  
  fDetDistCmd = new G4UIcmdWithADoubleAndUnit("/XFCT/det/distance",this);
  fDetDistCmd->SetGuidance("Set sample to detector window distance");
  fDetDistCmd->SetParameterName("DetDistance",false);
  fDetDistCmd->SetUnitCategory("Length");
  fDetDistCmd->AvailableForStates(G4State_Idle);
  
  
  fNanoConcCmd = new G4UIcmdWithADouble("/XFCT/det/setNanoConc",this);
  fNanoConcCmd->SetGuidance("Set concentration of nano particle soln");
  fNanoConcCmd->SetParameterName("NanoPartConc",false);
  fNanoConcCmd->AvailableForStates(G4State_Idle);
  
  fCollimRadCmd = new G4UIcmdWithADoubleAndUnit("/XFCT/det/setCollRadius",this);
  fCollimRadCmd->SetGuidance("Set radius of collimator");
  fCollimRadCmd->SetParameterName("CollimatorRadius",false);
  fCollimRadCmd->SetUnitCategory("Length");
  fCollimRadCmd->AvailableForStates(G4State_Idle);
  
  fCollimThickCmd = new G4UIcmdWithADoubleAndUnit("/XFCT/det/setCollThickness",this);
  fCollimThickCmd->SetGuidance("Set thickness of collimator");
  fCollimThickCmd->SetParameterName("CollimatorThickness",false);
  fCollimThickCmd->SetUnitCategory("Length");
  fCollimThickCmd->AvailableForStates(G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

XFCTDetectorMessenger::~XFCTDetectorMessenger()
{
  delete fFilter1MatCmd;
  delete fFilter2MatCmd;
  delete fThickFil1Cmd;
  delete fThickFil2Cmd;
  
  delete fFilterToSmplCmd;
  delete fDetDistCmd;
  
  delete fNanoConcCmd;
  
  delete fCollimRadCmd;
  delete fCollimThickCmd;
  
  delete fXFCTDirectory;
  delete fDetDirectory;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void XFCTDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if( command == fFilter1MatCmd )
   { fDetectorConstruction->SetFilter1Material(newValue);}

  if( command == fFilter2MatCmd )
   { fDetectorConstruction->SetFilter2Material(newValue);}
	   
  if( command == fThickFil1Cmd )
   { fDetectorConstruction->
	   SetFilter1Thickness(fThickFil1Cmd->GetNewDoubleValue(newValue));}
	   
  if( command == fThickFil2Cmd )
   { fDetectorConstruction->
	   SetFilter2Thickness(fThickFil2Cmd->GetNewDoubleValue(newValue));}
	   
  if( command == fFilterToSmplCmd )
   { fDetectorConstruction->
	   SetFilterToSampleDistance(fDetDistCmd->GetNewDoubleValue(newValue));}
	   
  if( command == fDetDistCmd )
   { fDetectorConstruction->
	   SetDetToSampleDistance(fDetDistCmd->GetNewDoubleValue(newValue));}
	   
  if( command == fNanoConcCmd )
   { fDetectorConstruction->
	   SetNanoParticleConc(fNanoConcCmd->GetNewDoubleValue(newValue));}
	   
  if( command == fCollimRadCmd )
   { fDetectorConstruction->
	   SetCollimRadius(fCollimRadCmd->GetNewDoubleValue(newValue));}

  if( command == fCollimThickCmd )
   { fDetectorConstruction->
	   SetCollimThickness(fCollimThickCmd->GetNewDoubleValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
