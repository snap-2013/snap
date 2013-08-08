#include "stdafx.h"

int main(int argc, char* argv[]) {
  Env = TEnv(argc, argv, TNotify::StdNotify);
  Env.PrepArgs(TStr::Fmt("GraphInfo. build: %s, %s. Time: %s", __TIME__, __DATE__, TExeTm::GetCurTm()));
  TExeTm ExeTm;
  Try
  const TStr InFNm = Env.GetIfArgPrefixStr("-i:", "as20graph.txt", "Input graph (one edge per line, tab/space separated)");
  const TBool OutputFormat = Env.GetIfArgPrefixBool("-d:", false, "Input Format (1 for Human-Friendly, 2 for Machine Friendly)");
  PUNGraph G = TSnap::LoadEdgeList<PUNGraph>(InFNm.CStr(),0,1);
  TGraphFeature Ori(G);
  Ori.SetAll();
  Ori.PrintVals(OutputFormat.Val);
  Catch
  printf("\nrun time: %s (%s)\n", ExeTm.GetTmStr(), TSecTm::GetCurTm().GetTmStr().CStr());
  return 0;
}
