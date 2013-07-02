#include "stdafx.h"

int ComputeKCore(const PUNGraph& G) {
  int cnt = 0;
  for(TUNGraph::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++)
    cnt = max(cnt, NI.GetOutDeg());
  THashSet <TInt> D[cnt+1];
  THash <TInt, TInt> deg;
  for(TUNGraph::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) { 
    TInt tmp = NI.GetOutDeg() - G->IsEdge(NI.GetId(), NI.GetId() );
    D[tmp.Val].AddKey(NI.GetId());
    deg.AddDat(NI.GetId()) = tmp;
  }
  int max_k = 0;
  for(int num_iters = 0;num_iters < G->GetNodes(); num_iters++)
    for(int i = 0; i < cnt; i++)
      if(D[i].Empty() == 0) {
        max_k = max(max_k, i);
        TInt a = *(D[i].BegI());
        D[i].DelKey(a);
        deg.AddDat(a.Val) = -1; // Hope overwriting works
        TUNGraph::TNodeI NI = G->GetNI(a.Val);
        for(int e = 0; e < NI.GetOutDeg(); e++) {
          TInt b = NI.GetOutNId(e);
          if(deg.GetDat(b) >= 0) {
            int Id = deg.GetKeyId(b); 
            D[deg[Id].Val].DelKey(b);
            deg[Id] = deg[Id] - 1;  //Hope the overwriting works
            D[deg[Id]].AddKey(b);
          }
        }
        break;
      }
  return max_k;
}

int main(int argc, char* argv[]) {
  PUNGraph G = TSnap::LoadEdgeList<PUNGraph>(argv[1],0,1);
  printf("%d\n",ComputeKCore(G));
  //Write your own function here
  return 0;
}
