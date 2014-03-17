// Compton Scattering and Pi0P/PiPN Photoproduction Event Generator
// Designed for use with the MAMI A2 Geant4 Simulation
//
// PiPN photoproduction event generator class
//
// Author - P. Martel
// Version - 23 October 2012

class PiPNGen : public BaseGen {
 protected:
  BasePart pPhoton, pTarget, pPiP, pRecoil;
  Float_t fPhotEk, fPhotEkCM;
  TLorentzVector lvMiss, lvPiP, lvPiPCM, lvReco, lvRecoCM;
  TH2F *hMiM;
  TH2F *hTpvtp, *hTpvtpi, *htpvtpi;
 public:
  PiPNGen(TString, TString, TString, Float_t, Float_t, Int_t, Int_t, Int_t);
  ~PiPNGen();
  void Init();
  Bool_t NewEvent(Float_t);
  void Reset();
  void SaveHists(TString);
};

PiPNGen::PiPNGen(TString sName, TString sTarget, TString sBase, Float_t fTargMass, Float_t fRecoMass, Int_t iRecoG3id, Int_t beamlo, Int_t beamhi) : BaseGen(sName, sTarget, sBase, beamlo, beamhi), pPhoton("Photon",0), pTarget("Target",fTargMass), pPiP("PiP",kMPI_MEV), pRecoil("Recoil",fRecoMass) {

  cout << "Constructing generator" << endl;

  // Setup Geant4 Ntuple

  Int_t ptag[2] = {iRecoG3id,8};

  InitNtuple(2,ptag);

  // Create some histograms for testing
  
  Int_t iMassCent = ((TMath::Nint(fRecoMass/50.))*50);

  hMiM = new TH2F("hMissM","Missing Mass vs Recoil K",100,0,500,60,iMassCent-150,iMassCent+150);
  hMiM->GetXaxis()->SetTitle("Recoil Kinetic Energy (MeV)");
  hMiM->GetYaxis()->SetTitle("Missing Mass (MeV)");

  hTpvtp = new TH2F("hTpvtp","#gamma p #rightarrow #pi^{+}n;#theta^{n} (deg);T^{n} (MeV)",180,0,180,180,0,180);
  hTpvtpi = new TH2F("hTpvtpi","#gamma p #rightarrow #pi^{+}n;#theta^{#pi^{+}}_{cm} (deg);T^{n} (MeV)",180,0,180,180,0,180);
  htpvtpi = new TH2F("htpvtpi","#gamma p #rightarrow #pi^{+}n;#theta^{#pi^{+}}_{cm} (deg);#theta^{n} (deg)",180,0,180,180,0,180);

};

PiPNGen::~PiPNGen(){

  cout << "Deleting generator" << endl;

  delete hMiM;
  delete hTpvtp;
  delete hTpvtpi;
  delete htpvtpi;

};

void PiPNGen::Init(){

  // Initialization for PiPN Photoproduction process

  if(bIsotW){
    cout << endl;
    cout << "--------------------------------------------------" << endl << endl;
    cout << "Using isotropic PiPN distributions" << endl;
  }
  
  else{
    cout << endl;
    cout << "--------------------------------------------------" << endl << endl;
    cout << "Loading PiPN Photoproduction data files" << endl;
    
    InitBase(sBaseName);

    cout << "E (MeV)\t\tMax CS (ub)\tTot CS (ub)" << endl;

    CrossGen();
  }

  cout << endl;
  cout << "--------------------------------------------------" << endl << endl;
  cout << "Running" << endl << endl;

  t1->Branch("Phot",&fPhotEk);
  t1->Branch("PhotCM",&fPhotEkCM);
  t1->Branch("PiP",&lvPiP);
  t1->Branch("PiPCM",&lvPiPCM);
  t1->Branch("Reco",&lvReco);
  t1->Branch("RecoCM",&lvRecoCM);

};

Bool_t PiPNGen::NewEvent(Float_t fBeamE){

  // Construct new event

  Reset();
  NewVertex();

  // Set initial state particles

  pPhoton.SetP4Lab(fBeamE,fBeamE,0,0);
  pTarget.SetP4Lab(pTarget.Mass,0,0,0);

  // Introduce a collision to determine energy, momentum, and direction
  // of final states particles

  Collision2B(pPhoton, pTarget, pPiP, pRecoil);

  // Check whether to reject the event based on the selected weighting

  if(Reject(fBeamE,pPiP.ThetaCM,pPiP.PhiCM)) return kFALSE;

  // For an incoherent process, use the previously determined directions of the
  // final state particles but determine the proper energies and momenta
  // through a spectral model
  
  if(bIncoh){
    SpecModel(pTarget);
    Collision2B(pPhoton, pTarget, pPiP, pRecoil);    
  }

  // Make sure the event is above threshold

  if(pPiP.Ener<pPiP.Mass){
    return kFALSE;
  }

  // Fill the particle histograms and construct the Geant4 Ntuple

  pPhoton.HistLab();
  pTarget.HistLab();
  pPiP.HistLab();
  pRecoil.HistLab();

  pPhoton.FillNtuple(var,3);
  pRecoil.FillNtuple(var,8);
  pPiP.FillNtuple(var,13);
  
  // Fill ntuple

  h1->Fill(var);
  
  // Fill tree

  fPhotEk = pPhoton.KEner;
  fPhotEkCM = pPhoton.KEnerCM;
  lvPiP = pPiP.P4;
  lvPiPCM = pPiP.P4CM;
  lvReco = pRecoil.P4;
  lvRecoCM = pRecoil.P4CM;

  t1->Fill();

  // Fill other test histograms

  lvMiss = (pPhoton.P4+pTarget.P4-pPiP.P4);

  hMiM->Fill(pRecoil.KEner,lvMiss.M());

  hTpvtp->Fill(pRecoil.Theta,pRecoil.KEner);
  hTpvtpi->Fill(pPiP.ThetaCM,pRecoil.KEner);
  htpvtpi->Fill(pPiP.ThetaCM,pRecoil.Theta);

  return kTRUE;
  
};

void PiPNGen::Reset(){

  // Reset particle information

  pPhoton.Reset();
  pTarget.Reset();
  pPiP.Reset();
  pRecoil.Reset();
  
};

void PiPNGen::SaveHists(TString sFile){

  // Write out test histograms

  cout << "Saving histograms" << endl;

  //Int_t i;

  TFile f1(sFile, "RECREATE", "MC_Hists_File");

  pPhoton.WriteHists();
  pTarget.WriteHists();
  pPiP.WriteHists();
  pRecoil.WriteHists();

  hMiM->Write();

  hTpvtp->Write();
  hTpvtpi->Write();
  htpvtpi->Write();

  f1.Close();

};
