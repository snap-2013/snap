void TGraphFeature::GrabStats(TFltV& V, TTuple <TFlt, EnumStatLast> & Stats) {
  Stats[Avg] = Stats[Std] = 0.0;
  V.Sort();
  Stats[Mn] = V[0];
  Stats[Md] = V[V.Len()/2];
  Stats[Mx] = V[V.Len()-1];
  for (int i = 0; i < V.Len(); i++) {
    Stats[Avg] += V[i];
    Stats[Std] += V[i] * V[i];
  }
  Stats[Avg] = Stats[Avg] / (double)V.Len();
  Stats[Std] = TMath::Sqrt(Stats[Std] / (double)V.Len() - Stats[Avg] * Stats[Avg]);
}

void TGraphFeature::SetDegPowerFit() {
  TVec<TFltPr> F(NumNodes);
  TVec <TInt> V(NumNodes);
  for (TUNGraph::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) { V[NI.GetOutDeg()]++; }
  for (int i = 0;i < NumNodes; i++) {
    if(V[i] > 0) { F.Add(TFltPr(i,V[i].Val)); }
  }
  TSpecFunc::PowerFit(F, A.Val, B.Val, SigA.Val, SigB.Val, Chi2.Val, R2.Val);
}

void TGraphFeature::SetDegCentr() {
  TVec <TFlt> V(NumNodes);
  for (TUNGraph::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
    V[NI.GetId()] = TSnap::GetDegreeCentr(Graph, NI.GetId());
  }
  GrabStats(V, DegC);
}

void TGraphFeature::SetBetwCentr() {  //If sampling is done, the sampled values are rescaled to reflect more accurately the actual Betweeness values
  double frac = min(1.0, (double)THRESHOLD_VAL / NumNodes);
  TIntFltH table;
  TSnap::GetBetweennessCentr(Graph,table,frac);
  TVec <TFlt> V(NumNodes);
  for (TIntFltH::TIter HI = table.BegI(); HI < table.EndI(); HI++) { V[HI.GetKey()] = HI.GetDat() / frac; }
  GrabStats(V, BetwC);
}

void TGraphFeature::SetClsCentr() {
  TVec <TFlt> V(NumNodes);
  for (TUNGraph::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
    V[NI.GetId()] = TSnap::GetClosenessCentr(Graph,NI.GetId());
  }
  GrabStats(V, ClsC);
}

void TGraphFeature::SetKCore() {
  TKCore <PUNGraph> KG (Graph);
  KCoreSize = 0;
  while(1) {
    KG.GetNextCore();
    if(KG.GetCoreNodes() == 0) break;
    KCoreSize = KCoreSize + 1;
  }
  FracKCore = (double)(TSnap::GetKCore(Graph, KCoreSize))->GetNodes() / (double)NumNodes;
}

void TGraphFeature::SetCCF() {
  TVec <TFlt> V(NumNodes);
  for (TUNGraph::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
    V[NI.GetId()] = TSnap::GetNodeClustCf(Graph,NI.GetId());
  }
  GrabStats(V, CCF);
  TFltPrV DegCCfV;
  CCF_Global = TSnap::GetClustCf(Graph, DegCCfV, ClosedTriads, OpenTriads);
  FracTriads = (double)ClosedTriads / (double) (OpenTriads + ClosedTriads);
}

void TGraphFeature::SetDiam() {
  TSnap::GetBfsEffDiam(Graph, min(THRESHOLD_VAL, NumNodes.Val), 0, EffDiam.Val, FullDiam.Val); //REMEMBER TO SET THRESHOLD VAL!!!
}

void TGraphFeature::SetFracLCC() {
  TIntPrV V;
  TSnap::GetWccSzCnt(Graph, V);
  TIntPr tmp = V[V.Len()-1];
  FracLCC[First] = (double)tmp.Val1 / NumNodes;
  if (tmp.Val2 > 1) { FracLCC[Second] = FracLCC[First]; }
  else if(V.Len() == 1) FracLCC[Second] = 0.0;
  else { FracLCC[Second] = (double)V[V.Len()-2].Val1 / NumNodes; }
}

void TGraphFeature::SetEgnVal() {
  TFltV V;
  TSnap::GetEigVals(Graph, 2, V);
  EgnVal[First] = V[0];
  EgnVal[Second] = V[1];
}

void TGraphFeature::SetAssty() {  
  AsstyCor = TSnap::GetAsstyCor(Graph);
}

void TGraphFeature::SetAll() {
  SetDegPowerFit();
  SetDegCentr();
  SetBetwCentr();
  SetClsCentr();
  SetKCore();
  SetCCF();
  SetDiam();
  SetFracLCC();
  SetEgnVal();
  SetAssty();
}
