#include "stdafx.h"

int main (int argc, char *argv[]) {
  PUNGraph G = TSnap::LoadEdgeList<PUNGraph>(argv[1],0,1);
  TSnap::RenumberNodes(G);
  TSnap::RemoveSelfLoops(G);
  TVec <int64> A, B;
  TSnap::GetMotifCount(G,3,A,G->GetNodes()); //Replace G->GetNodes() with a smaller integer value if one wants to do sampling
  TSnap::GetMotifCount(G,4,B,G->GetEdges()); //Replace G->GetEdges() with a smaller integer value if one wants to do sampling
  printf("%s,%d,%d",argv[1],G->GetNodes(),G->GetEdges());
  for(int i=0;i<2;++i) printf(",%lld",A[i]);
  for(int i=0;i<6;++i) printf(",%lld",B[i]);
  printf("\n");
}

