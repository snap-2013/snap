namespace TSnap {

void GetMotifCount(const PUNGraph& G, const int MotifSize, TVec <int64> & MotifV, const int num) {
  if (MotifSize == 3) {
    MotifV = TVec <int64> (2);
    MotifV.PutAll(0);
    TSnap::GetTriads(G,MotifV[mtThreeClosed],MotifV[mtThreeOpen],num);
  }
  else {
    MotifV = TVec <int64> (6);
    MotifV.PutAll(0);
    TIntPrV V(G->GetEdges(), 0);
    for (TUNGraph::TEdgeI EI = G->BegEI(); EI < G->EndEI(); EI++) {
      V.Add(TIntPr(EI.GetSrcNId(), EI.GetDstNId()));
    }
    TRnd blargh;
    V.Shuffle(blargh);
    for (int z = 0; z < num; z++) {
      int SrcNId = V[z].Val1.Val, DstNId = V[z].Val2.Val;
      TUNGraph::TNodeI SrcNI = G->GetNI(SrcNId), DstNI = G->GetNI(DstNId);
      TIntV SrcV(SrcNI.GetOutDeg(),0), DstV(DstNI.GetOutDeg(),0), BothV(min(SrcNI.GetOutDeg(), DstNI.GetOutDeg()),0);
      SrcV.Clr(0,-1);
      DstV.Clr(0,-1);
      BothV.Clr(0,-1);
      //Grouping the vertices into sets
      for (int e = 0; e < SrcNI.GetOutDeg(); e++) {
        if (G->IsEdge(DstNId, SrcNI.GetOutNId(e)) ) { BothV.Add(SrcNI.GetOutNId(e)); }
        else { SrcV.Add(SrcNI.GetOutNId(e)); }
      }
      for (int e = 0; e < DstNI.GetOutDeg(); e++) {
        if (G->IsEdge(SrcNId, DstNI.GetOutNId(e)) == 0) { DstV.Add(DstNI.GetOutNId(e)); }
      }
      //Compute Motif 0 and 1
      for (int i = 0; i < SrcV.Len(); i++) {
        for (int j = 0; j < DstV.Len(); j++) {
          if (G->IsEdge(SrcV[i], DstV[j]) ) { MotifV[mfFourSquare]++; }
          else MotifV[mfFourLine]++;
        }
      }
      //Compute Motif 2 and 3
      for (int i = 0; i < SrcV.Len(); i++) {
        for (int j = i + 1; j < SrcV.Len(); j++) {
          if (G->IsEdge(SrcV[i], SrcV[j]) ) { MotifV[mfFourTriangleEdge]++; }
          else MotifV[mfFourStar]++;
        }
      }
      for (int i = 0; i < DstV.Len(); i++) {
        for (int j = i + 1; j < DstV.Len(); j++) {
          if (G->IsEdge(DstV[i], DstV[j]) ) { MotifV[mfFourSquare]++; }
          else MotifV[mfFourStar]++;
        }
      }
      //Compute Motif 4 and 5
      for (int i = 0; i < BothV.Len(); i++) {
        for (int j = i + 1; j < BothV.Len(); j++) {
          if (G->IsEdge(BothV[i], BothV[j]) ) { MotifV[mfFourComplete]++; }
          else MotifV[mfFourSquareDiag]++;
        }
      }
    }
    MotifV[mfFourSquare] /= 4ll;
    MotifV[mfFourStar] /= 3ll;
    MotifV[mfFourSquareDiag] /= 6ll;
  }
}

}; //namespace TSnap
