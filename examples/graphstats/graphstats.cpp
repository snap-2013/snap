#include "stdafx.h"

int main(int argc, char* argv[]) {
  PUNGraph G = TSnap::LoadEdgeList<PUNGraph>(argv[1],0,1);
  TGraphFeature Ori(G);
  Ori.SetAll();
  printf("Default Output Format: \n");
  Ori.PrintVals();
  printf("Machine-Friendly Output Format: \n");
  Ori.PrintVals(1);
  return 0;
}
