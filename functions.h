void SpecModel(TLorentzVector& pki){

  Float_t kpi, kth, kph, prob, func;
  Float_t N_s = 0.52;
  Float_t N_p = 1.65;
  Float_t k_max = 300.0;
  Float_t A[2] = {43.5, 13500.0};
  Float_t coeff[2] = {0.0, 2.0};
  Float_t sigma[2] = {90.0, 75.0};
  Float_t Emin[2] = {26.0, 16.0};
  Float_t Emax[2] = {50.0, 26.0};
  Float_t sf_max[2] = {0.26, 0.93};

  Int_t L;

  if((gRandom->Rndm()) < N_s/(N_s+N_p)) L = 0;
  else L = 1;
  
  prob = 1;
  func = 0;

  while(prob>func){
    kpi = k_max*gRandom->Rndm();
    prob = sf_max[L]*gRandom->Rndm();
    func = A[L]*pow((kpi/1000.0),(coeff[L]+2))*exp(-Sqr(kpi/sigma[L])/2);
  }

  kth = acos(-1+2*gRandom->Rndm());
  kph = kPI*(-1+2*gRandom->Rndm());
    
  pki.SetPxPyPzE(1,1,1,kMP_MEV-(Emin[L]+(Emax[L]-Emin[L])*gRandom->Rndm()));
  pki.SetRho(kpi);
  pki.SetTheta(kth);
  pki.SetPhi(kph);

}

TString GenNames( Int_t npart, Int_t* ptag) {
  
  Int_t i, j;
  
  TString pstr[] = {"Px", "Py", "Pz", "Pt", "En"};
  TString beam = "X_vtx:Y_vtx:Z_vtx:Px_bm:Py_bm:Pz_bm:Pt_bm:En_bm";
  TString particles;
  TString names;
  
  for ( i = 0; i < npart; i++) {
    for ( j = 0; j < 5; j++) {
      particles.Append( pstr[j]);
      if ( ( i == (npart-1)) && ( j == 4))
	particles.Append( Form( "_l%02d%02d", i+1, ptag[i]));
      else
	particles.Append( Form( "_l%02d%02d:", i+1, ptag[i]));
    }
  }
  
  names = beam + ":" + particles;
  
  return( names);
}

Float_t Sqr(Float_t x)
{
  return(x*x);
}
