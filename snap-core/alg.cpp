namespace TSnap {

double GetAsstyCor(const PUNGraph& Graph) {
  TIntFltH deg(Graph->GetNodes()), deg_sq(Graph->GetNodes());
  for (TUNGraph::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
    deg.AddDat(NI.GetId()) = NI.GetOutDeg(); 
    deg_sq.AddDat(NI.GetId()) = NI.GetOutDeg() * NI.GetOutDeg();
  }
  double m = Graph->GetEdges(), num1 = 0.0, num2 = 0.0, den1 = 0.0;
  for (TUNGraph::TEdgeI EI = Graph->BegEI(); EI < Graph->EndEI(); EI++) {
    double t1 = deg.GetDat(EI.GetSrcNId()).Val, t2 = deg.GetDat(EI.GetDstNId()).Val;
    num1 += t1 * t2;
    num2 += t1 + t2;
    den1 += deg_sq.GetDat(EI.GetSrcNId()).Val + deg_sq.GetDat(EI.GetDstNId()).Val;
  }
  num1 /= m;
  den1 /= (2.0 * m);
  num2 = (num2 / (2.0 * m)) * (num2 / (2.0 * m));
  return (num1 - num2) / (den1 - num2);
}

}; //namespace TSnap
