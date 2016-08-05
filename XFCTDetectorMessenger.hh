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
// $Id: B2aDetectorMessenger.hh 69706 2013-05-13 09:12:40Z gcosmo $
//
/// \file B2aDetectorMessenger.hh
/// \brief Definition of the B2aDetectorMessenger class

#ifndef XFCTDetectorMessenger_h
#define XFCTDetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class XFCTDetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithADouble;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// Messenger class that defines commands for B2aDetectorConstruction.
///
/// It implements commands:
/// - /XFCT/det/setFilter1Material name
/// - /XFCT/det/setFilter2Material name
/// - /XFCT/det/distance value unit

class XFCTDetectorMessenger: public G4UImessenger
{
  public:
    XFCTDetectorMessenger(XFCTDetectorConstruction* );
    virtual ~XFCTDetectorMessenger();
    
    virtual void SetNewValue(G4UIcommand*, G4String);
    
  private:
    XFCTDetectorConstruction*  fDetectorConstruction;

    G4UIdirectory*           		fXFCTDirectory;
    G4UIdirectory*           		fDetDirectory;

    G4UIcmdWithAString*      		fFilter1MatCmd;
    G4UIcmdWithAString*      		fFilter2MatCmd;

    
    G4UIcmdWithADoubleAndUnit*		fThickFil1Cmd;
    G4UIcmdWithADoubleAndUnit*		fThickFil2Cmd;
    G4UIcmdWithADoubleAndUnit*		fFilterToSmplCmd;
    G4UIcmdWithADoubleAndUnit*		fDetDistCmd;
    G4UIcmdWithADoubleAndUnit*		fCollimRadCmd;
    G4UIcmdWithADoubleAndUnit*		fCollimThickCmd;
    
    G4UIcmdWithADouble*				fNanoConcCmd;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
