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
// $Id: XFCTPrimaryGeneratorAction.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file XFCTPrimaryGeneratorAction.cc
/// \brief Implementation of the XFCTPrimaryGeneratorAction class

#include "XFCTPrimaryGeneratorAction.hh"
#include "XFCTPrimaryGeneratorMessenger.hh"
#include "XFCTRunAction.hh"


#include "G4RunManager.hh"
#include "G4MTRunManager.hh"

#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

#include "G4PhysicalConstants.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "Randomize.hh"

//----------------------------------------------------------------------
//		Default particle kinetics: Point source with zero divergence
//----------------------------------------------------------------------
XFCTPrimaryGeneratorAction::XFCTPrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  coneBeam("on")
{
  runAction = 0;
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);

  //create a messenger for this class
  gunMessenger = new XFCTPrimaryGeneratorMessenger(this);

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
    = particleTable->FindParticle(particleName="gamma");
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  fParticleGun->SetParticleEnergy(125.*keV);
  
  G4double x0, y0, z0;
  x0 = y0 = 0.*cm;
  z0 = -3.55*cm;
  
  fParticleGun->SetParticlePosition(G4ThreeVector(x0, y0, z0));

  G4cout << "XFCTPrimaryGeneratorAction created" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

XFCTPrimaryGeneratorAction::~XFCTPrimaryGeneratorAction()
{
  delete fParticleGun;
  delete gunMessenger;
}

//----------------------------------------------------------------------
//				Real particle kinetics: Cone beam source
//----------------------------------------------------------------------


void XFCTPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //retrieve runAction, if not done
  if (!runAction)
    {
      //Sequential runaction
      if (G4RunManager::GetRunManager()->GetRunManagerType() == 
	  G4RunManager::sequentialRM)
	runAction = static_cast<const XFCTRunAction*>
	  (G4RunManager::GetRunManager()->GetUserRunAction());  
      else //MT master runaction
	runAction = static_cast<const XFCTRunAction*>
	  (G4MTRunManager::GetMasterRunManager()->GetUserRunAction());  
      if (!runAction)
	G4cout << "Something wrong here!" << G4endl;
    }
 
  //this function is called at the begining of event 
  
  if (coneBeam == "on")
    {
	  //G4cout <<" Cone beam is used to generate particle" << G4endl;
	  G4double z0 = -3.55*cm;
      G4double focusRadius = 0.55*0.5*cm;
      
      //theta in [0;pi/2]
      G4double theta = (pi/9)*G4UniformRand();
      
      //phi in [-pi;pi]
      G4double phi = (G4UniformRand()*2*pi)- pi;
      G4double x = focusRadius*std::sin(theta)*std::sin(phi);
      G4double y = focusRadius*std::sin(theta)*std::cos(phi);
      G4double z = -(focusRadius*std::cos(theta));
      
      fParticleGun->SetParticlePosition(G4ThreeVector(x,y,z0));
      
      G4double Xdim = focusRadius;
      G4double Ydim = focusRadius;
      
      G4double Dx = Xdim*(G4UniformRand()-0.5);
      
      G4double Dy = Ydim*(G4UniformRand()-0.5);
      
      fParticleGun->SetParticleMomentumDirection(G4ThreeVector(-x+Dx,-y+Dy,-z));
      
    }
 
  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....









//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

