void TGraphFeature::GrabStats(TFltV& V, TTuple <TFlt, fsEnumStatLast> & Stats) {
  Stats[fsAvg] = Stats[fsStd] = 0.0;
  V.Sort();
  Stats[fsMn] = V[0];
  Stats[fsMd] = V[V.Len()/2];
  Stats[fsMx] = V[V.Len()-1];
  for (int i = 0; i < V.Len(); i++) {
    Stats[fsAvg] += V[i];
    Stats[fsStd] += V[i] * V[i];
  }
  Stats[fsAvg] = Stats[fsAvg] / (double)V.Len();
  Stats[fsStd] = TMath::Sqrt(Stats[fsStd] / (double)V.Len() - Stats[fsAvg] * Stats[fsAvg]);
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
  double frac = min(1.0, (double)MaxSamplingSize / NumNodes);
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
  TSnap::GetBfsEffDiam(Graph, min(MaxSamplingSize, NumNodes.Val), 0, EffDiam.Val, FullDiam.Val);
}

void TGraphFeature::SetFracLCC() {
  TIntPrV V;
  TSnap::GetWccSzCnt(Graph, V);
  TIntPr tmp = V[V.Len()-1];
  FracLCC[fsFirst] = (double)tmp.Val1 / NumNodes;
  if (tmp.Val2 > 1) { FracLCC[fsSecond] = FracLCC[fsFirst]; }
  else if(V.Len() == 1) FracLCC[fsSecond] = 0.0;
  else { FracLCC[fsSecond] = (double)V[V.Len()-2].Val1 / NumNodes; }
}

void TGraphFeature::SetEgnVal() {
  TFltV V;
  TSnap::GetEigVals(Graph, 2, V);
  EgnVal[fsFirst] = V[0];
  EgnVal[fsSecond] = V[1];
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
