void TGraphFeature::SetAvgStd(TFltV& V, TFlt & avg, TFlt & std) {
  avg = std = 0.0;
  for (int i = 0; i < V.Len(); i++) {
    avg += V[i];
    std += V[i] * V[i];
  }
  avg = avg / (double)V.Len();
  std = TMath::Sqrt(std / (double)V.Len() - avg * avg);
}

void TGraphFeature::SetFracDeg() {
  TInt t1 = 0, t2 = 0;
  TFlt AvgDeg = 0;
  for (TUNGraph::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) { AvgDeg += NI.GetOutDeg(); }
  AvgDeg = AvgDeg / (double)NumNodes;
  for (TUNGraph::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
    if (NI.GetOutDeg() >= AvgDeg.Val) t1 = t1 + 1;
    if (NI.GetOutDeg() >= 2.0 * AvgDeg.Val) t2 = t2 + 1;
  }
  StatV[fsFracDeg] = (double)t1 / (double)NumNodes;
  StatV[fsFrac2Deg] = (double)t2 / (double)NumNodes; 
}

void TGraphFeature::SetBiCon() {
  TIntPrV V;
  TSnap::GetBiConSzCnt(Graph, V);
  StatV[fsFracBiCon] = 0.0;
  StatV[fsFracBiCon] = (double)V[V.Len()-1].Val1.Val / (double)NumNodes; 
}

void TGraphFeature::SetDegCentr() {
  TVec <TFlt> V(NumNodes);
  for (TUNGraph::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
    V[NI.GetId()] = TSnap::GetDegreeCentr(Graph, NI.GetId());
  }
  SetAvgStd(V, StatV[fsAvgDegC], StatV[fsStdDegC]);
}

void TGraphFeature::SetBetwCentr() {  //If sampling is done, the sampled values are rescaled to reflect more accurately the actual Betweeness values
  double frac = min(1.0, (double)MaxSamplingSize / NumNodes);
  TIntFltH table;
  TSnap::GetBetweennessCentr(Graph,table,frac);
  TVec <TFlt> V(NumNodes);
  for (TIntFltH::TIter HI = table.BegI(); HI < table.EndI(); HI++) { V[HI.GetKey()] = HI.GetDat() / frac; }
  TFlt minVal = V[0].Val, maxVal = V[0].Val;
  for (int i = 0; i < V.Len(); i++) {
    minVal = min(minVal.Val, V[i].Val);
    maxVal = max(maxVal.Val, V[i].Val);
  }
  for (int i = 0; i < V.Len(); i++) { V[i] = (V[i].Val - minVal.Val) / (maxVal.Val - minVal.Val); }
  SetAvgStd(V, StatV[fsAvgBetwC], StatV[fsStdBetwC]);
}

void TGraphFeature::SetClsCentr() {
  int num = min(MaxSamplingSize, NumNodes.Val);
  TVec <TFlt> V(num);
  for (int i = 0; i < num; i++) {
    V[i] = TSnap::GetClosenessCentr(Graph, Graph->GetRndNId());
  }
  SetAvgStd(V, StatV[fsAvgClsC], StatV[fsStdClsC]);
}

void TGraphFeature::SetKCore() {
  TKCore <PUNGraph> KG (Graph);
  StatV[fsFrac2Core] = StatV[fsFrac3Core] = StatV[fsFrac4Core] = 0.0;
  TInt tmp = 0;
  while(1) {
    KG.GetNextCore();
    if(KG.GetCoreNodes() == 0) break;
    tmp = tmp + 1;
    if(tmp.Val == 2) { StatV[fsFrac2Core] = (double)KG.GetCoreNodes() / (double)NumNodes; }
    else if(tmp.Val == 3) { StatV[fsFrac3Core] = (double)KG.GetCoreNodes() / (double)NumNodes; }
    else if(tmp.Val == 4) { StatV[fsFrac4Core] = (double)KG.GetCoreNodes() / (double)NumNodes; }
  }
  StatV[fsKCoreSize] = tmp;
  StatV[fsFracKCore] = (double)(TSnap::GetKCore(Graph, tmp))->GetNodes() / (double)NumNodes;
}

void TGraphFeature::SetCCF() {
  TVec <TFlt> V(NumNodes);
  for (TUNGraph::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
    V[NI.GetId()] = TSnap::GetNodeClustCf(Graph,NI.GetId());
  }
  SetAvgStd(V, StatV[fsAvgCCF], StatV[fsStdCCF]);
}

void TGraphFeature::SetDiam() {
  int tmp;
  TSnap::GetBfsEffDiam(Graph, min(MaxSamplingSize, NumNodes.Val), 0, StatV[fsEffDiam].Val, tmp);
  StatV[fsFullDiam] = tmp;
}

void TGraphFeature::SetFracLCC() {
  TIntPrV V;
  TSnap::GetWccSzCnt(Graph, V);
  TIntPr tmp = V[V.Len()-1];
  StatV[fsFracLCC1] = (double)tmp.Val1 / (double)NumNodes;
  if (tmp.Val2 > 1) { StatV[fsFracLCC2] = StatV[fsFracLCC1]; }
  else if(V.Len() == 1) StatV[fsFracLCC2] = 0.0;
  else { StatV[fsFracLCC2] = (double)V[V.Len()-2].Val1 / (double)NumNodes; }
}

void TGraphFeature::SetEgnVal() {
  TFltV V;
  TSnap::GetEigVals(Graph, 2, V);
  StatV[fsEgnVal1] = V[0];
  StatV[fsEgnVal2] = V[1];
  if(StatV[fsEgnVal2].Val != StatV[fsEgnVal2].Val) { StatV[fsEgnVal2] = 0.0; } // Check for NaN
}

void TGraphFeature::SetAssty() {  
  StatV[fsAsstyCor] = TSnap::GetAsstyCor(Graph);
}

void TGraphFeature::SetMotif() {
  TVec <int64> A, B;
  TInt num = min(NumEdges.Val, MaxMotSampleSize);
  TSnap::GetMotifCount(Graph, 3, A, min(NumNodes.Val, MaxMotSampleSize));
  TSnap::GetMotifCount(Graph, 4, B, min(NumEdges.Val, MaxMotSampleSize));
  int64 t1 = 0, t2 = 0;
  for (int i = 0; i < 2; i++) t1 += A[i];
  MotV[fm3closed] = (t1 == 0)?0.0:((double)A[0] / t1);
  MotV[fm3open] = (t1 == 0)?0.0:((double)A[1] / t1);
  for (int i = 0; i < 6; i++) t2 += B[i];
  for (int it = fm4line; it < fmLast; it++) {
    MotV[it] = (t2 == 0)?0.0:((double)B[it - fm4line] / t2);
  }
}

void TGraphFeature::SetAll() {
  SetFracDeg();
  SetBiCon();
  SetDegCentr();
  SetBetwCentr();
  SetClsCentr();
  SetKCore();
  SetCCF();
  SetDiam();
  SetFracLCC();
  SetEgnVal();
  SetAssty();
  SetMotif();
}



