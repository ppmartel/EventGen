// Compton Scattering and Pi0 Photoproduction Event Generator
// Designed for use with the MAMI A2 Geant4 Simulation
//
// Base particle class for use in the event generator
//
// Author - P. Martel
// Version - 23 June 2011

class BasePart {
 protected:
  TString sPartName;
  void FillCM();
  void FillLab();
 public:
  BasePart(TString, Float_t);
  ~BasePart();
  TLorentzVector P4, P4CM;
  Float_t Mass, Ener, KEner, Mom, Theta, Phi;
  Float_t EnerCM, KEnerCM, MomCM, ThetaCM, PhiCM;
  TH1F *hKEner, *hMom, *hTheta, *hPhi, *hMass;
  void BoostCM(TVector3);
  void BoostLab(TVector3);
  void HistCM();
  void HistLab();
  void Reset();
  void SetP4CM(Float_t, Float_t, Float_t, Float_t);
  void SetP4Lab(Float_t, Float_t, Float_t, Float_t);
  void SetP4CM(Float_t, Float_t);
  void SetP4Lab(Float_t, Float_t);
  void RotateCM(Float_t);
  void RotateLab(Float_t);
  void RotatePhi(Float_t);
  TString WhichDet();
  void WriteHists();
  void FillNtuple(Float_t*, Int_t);
};

BasePart::BasePart(TString name, Float_t mass){

  sPartName = name;
  Mass = mass;

  cout << ("Constructing "+sPartName) << endl;

  // Create some histograms for testing
  
  hKEner = new TH1F(sPartName+"KEner",sPartName+" Kinetic E (MeV)",100,0,500);
  hKEner->GetXaxis()->SetTitle("Kinetic Energy (MeV)");
  hKEner->GetYaxis()->SetTitle("# of Events");

  hMom = new TH1F(sPartName+"Mom",sPartName+" Momentum (MeV)",100,0,500);
  hMom->GetXaxis()->SetTitle("Momentum (MeV)");
  hMom->GetYaxis()->SetTitle("# of Events");

  hTheta = new TH1F(sPartName+"Theta",sPartName+" Theta (deg)",180,0,180);
  hTheta->GetXaxis()->SetTitle("Theta (deg)");
  hTheta->GetYaxis()->SetTitle("# of Events");

  hPhi = new TH1F(sPartName+"Phi",sPartName+" Phi (deg)",360,-180,180);
  hPhi->GetXaxis()->SetTitle("Phi (deg)");
  hPhi->GetYaxis()->SetTitle("# of Events");

  Int_t MassCent = ((TMath::Nint(Mass/50.))*50);

  hMass = new TH1F(sPartName+"Mass",sPartName+" Mass (MeV)",60,MassCent-150,MassCent+150);
  hMass->GetXaxis()->SetTitle("Mass (MeV)");
  hMass->GetYaxis()->SetTitle("# of Events");

};

BasePart::~BasePart(){

  cout << ("Deleting "+sPartName) << endl;

  delete hKEner;
  delete hMom;
  delete hTheta;
  delete hPhi;
  delete hMass;

};

void BasePart::FillCM(){

  // Fill CM variables using the already set Lorentz vector

  EnerCM = P4CM.E();
  KEnerCM = (EnerCM-Mass);
  MomCM = P4CM.Rho();
  ThetaCM = (P4CM.Theta()*kR2D);
  PhiCM = (P4CM.Phi()*kR2D);

};

void BasePart::FillLab(){

  // Fill lab variables using the already set Lorentz vector

  Ener = P4.E();
  KEner = (Ener-Mass);
  Mom = P4.Rho();
  Theta = (P4.Theta()*kR2D);
  Phi = (P4.Phi()*kR2D);

};

void BasePart::BoostCM(TVector3 v){

  // Boost INTO the CM frame, and fill its variables

  P4CM = P4;
  P4CM.Boost(v);
  FillCM();

};

void BasePart::BoostLab(TVector3 v){

  // Boost INTO the lab frame, and fill its variables

  P4 = P4CM;
  P4.Boost(v);
  FillLab();

};

void BasePart::HistCM(){

  // Fill CM histograms

  hKEner->Fill(KEnerCM);
  hMom->Fill(MomCM);
  hTheta->Fill(ThetaCM);
  hPhi->Fill(PhiCM);
  hMass->Fill(Mass);

};

void BasePart::HistLab(){

  // Fill lab histograms

  hKEner->Fill(KEner);
  hMom->Fill(Mom);
  hTheta->Fill(Theta);
  hPhi->Fill(Phi);
  hMass->Fill(Mass);

};

void BasePart::Reset(){

  // Reset CM and lab Lorentz vectors and variables

  P4CM.SetPxPyPzE(0,0,0,0);
  P4.SetPxPyPzE(0,0,0,0);
  FillCM();
  FillLab();

};

void BasePart::SetP4CM(Float_t ener,Float_t mom,Float_t theta,Float_t phi){

  // Set the CM Lorentz Vector with the supplied values

  P4CM.SetPxPyPzE(1,1,1,ener);
  P4CM.SetRho(mom);
  P4CM.SetTheta(theta*kD2R);
  P4CM.SetPhi(phi*kD2R);
  FillCM();

};

void BasePart::SetP4Lab(Float_t ener,Float_t mom,Float_t theta,Float_t phi){

  // Set the lab Lorentz Vector with the supplied values

  P4.SetPxPyPzE(1,1,1,ener);
  P4.SetRho(mom);
  P4.SetTheta(theta*kD2R);
  P4.SetPhi(phi*kD2R);
  FillLab();

};

void BasePart::SetP4CM(Float_t ener, Float_t mom){

  // Set the CM Lorentz Vector with the supplied values
  // providing an isotropic distribution

  Float_t theta = acos(-1+2*gRandom->Rndm());
  Float_t phi = kPI*(-1+2*gRandom->Rndm());
  SetP4CM(ener, mom, (theta*kR2D), (phi*kR2D));

};

void BasePart::SetP4Lab(Float_t ener, Float_t mom){

  // Set the lab Lorentz Vector with the supplied values
  // providing an isotropic distribution

  Float_t theta = acos(-1+2*gRandom->Rndm());
  Float_t phi = kPI*(-1+2*gRandom->Rndm());
  SetP4Lab(ener, mom, (theta*kR2D), (phi*kR2D));

};

void BasePart::RotateCM(Float_t phi){

  // Rotate CM frame about z-axis by phi

  P4CM.RotateZ(kD2R*phi);
  FillCM();

};

void BasePart::RotateLab(Float_t phi){

  // Rotate lab frame about z-axis by phi

  P4.RotateZ(kD2R*phi);
  FillLab();

};

void BasePart::RotatePhi(Float_t phi){

  // Rotate both frames about z-axis by phi

  RotateCM(phi);
  RotateLab(phi);

};

TString BasePart::WhichDet(){

  // Simple determination of which detector the particle will hit

  TString sDet;
  if(Theta<=20) sDet = "TAPS";
  else if(Theta<=160) sDet = "CB";
  else sDet = "Out";

  return sDet;

};

void BasePart::WriteHists(){

  // Write out the particle histograms

  hKEner->Write();
  hMom->Write();
  hTheta->Write();
  hPhi->Write();
  hMass->Write();

};

void BasePart::FillNtuple(Float_t* var, Int_t index){

  // Construct the Geant4 Ntuple for this particle

  var[index++] = sin(Theta*kD2R)*cos(Phi*kD2R);
  var[index++] = sin(Theta*kD2R)*sin(Phi*kD2R);
  var[index++] = cos(Theta*kD2R);
  var[index++] = Mom/1000;
  var[index++] = Ener/1000;

};
