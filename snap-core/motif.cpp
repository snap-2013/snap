namespace TSnap {

/*
 * For MotifSize = 3, 0 is ClosedTriads, 1 is OpenTriads
 *
 * For MotifSize = 4, we use the following convention to label the motifs
 *
 * NOTE: The pair (a,b) represents an edge between a and b 
 *
 * 0 : (1,2), (2,3), (3,4)
 *
 * 1 : (1,2), (2,3), (3,4), (1,4)
 *
 * 2 : (1,2), (1,3), (1,4)
 *
 * 3 : (1,2), (1,3), (1,4), (2,3)
 *
 * 4 : (1,2), (1,3), (1,4), (2,3), (3,4)
 *
 * 5 : (1,2), (1,3), (1,4), (2,3), (3,4), (2,4)
 *
 */

void GetMotifCount(const PUNGraph& G, const int MotifSize, TVec <int64> & MotifV) {
  if (MotifSize == 3) {
    MotifV = TVec <int64> (2);
    MotifV.PutAll(0);
    TSnap::GetTriads(G,MotifV[0],MotifV[1],-1);
  }
  else {
    MotifV = TVec <int64> (6);
    MotifV.PutAll(0);
    int MxVSize = 0;
    for (TUNGraph::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
      MxVSize = max(MxVSize, NI.GetOutDeg());
    }
    TIntV SrcV(MxVSize), DstV(MxVSize), BothV(MxVSize);
    for (TUNGraph::TEdgeI EI = G->BegEI(); EI < G->EndEI(); EI++) {
      TUNGraph::TNodeI SrcNI = G->GetNI(EI.GetSrcNId()), DstNI = G->GetNI(EI.GetDstNId());
      int SrcNId = SrcNI.GetId(), DstNId = DstNI.GetId();
    //  TIntV SrcV(SrcNI.GetOutDeg(),0), DstV(DstNI.GetOutDeg(),0), BothV(min(SrcNI.GetOutDeg(), DstNI.GetOutDeg()),0);
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
          if (G->IsEdge(SrcV[i], DstV[j]) ) { MotifV[1]++; }
          else MotifV[0]++;
        }
      }
      //Compute Motif 2 and 3
      for (int i = 0; i < SrcV.Len(); i++) {
        for (int j = i + 1; j < SrcV.Len(); j++) {
          if (G->IsEdge(SrcV[i], SrcV[j]) ) { MotifV[3]++; }
          else MotifV[2]++;
        }
      }
      for (int i = 0; i < DstV.Len(); i++) {
        for (int j = i + 1; j < DstV.Len(); j++) {
          if (G->IsEdge(DstV[i], DstV[j]) ) { MotifV[3]++; }
          else MotifV[2]++;
        }
      }
      //Compute Motif 4 and 5
      for (int i = 0; i < BothV.Len(); i++) {
        for (int j = i + 1; j < BothV.Len(); j++) {
          if (G->IsEdge(BothV[i], BothV[j]) ) { MotifV[5]++; }
          else MotifV[4]++;
        }
      }
    }
    MotifV[1] /= 4ll;
    MotifV[2] /= 3ll;
    MotifV[5] /= 6ll;
  }
}

}; //namespace TSnap
