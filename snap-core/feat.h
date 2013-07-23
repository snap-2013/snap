const int MaxSamplingSize = 30000;
const int MaxMotSampleSize = 70000;

typedef enum {
  fsAvgCCF=0, fsAvgDegC, fsAvgBetwC, fsAvgClsC, fsFullDiam, fsEffDiam,
  fsFracLCC1, fsFracLCC2, fsFracDeg, fsFrac2Deg, fsFracBiCon, 
  fsKCoreSize, fsFracKCore, fsFrac2Core, fsFrac3Core, fsFrac4Core, 
  fsEgnVal1, fsEgnVal2, fsAsstyCor,  
  fsStdCCF, fsStdDegC, fsStdBetwC, fsStdClsC, 
  fsLast
}TGraphFeatureStats;

typedef enum {
  fm3closed=0, fm3open,
  fm4line, fm4cycle, fm4star, fm4triedge, fm4edgemissing, fm4complete,
  fmLast
}TGraphFeatureMotifs;

class TGraphFeature {
private:
  PUNGraph Graph;
public:
  TInt NumNodes, NumEdges;
  TFltV StatV, MotV;
private:
  void SetAvgStd (TFltV& V, TFlt & avg, TFlt & std);                 //Computes Average and Standard Deviation 
public:
  TGraphFeature (PUNGraph& _Graph) : Graph(_Graph) {                 //Constructor. Nodes relabelled from 0 to N-1 and self-loops removed.
    TSnap::RenumberNodes(Graph);                                     
    TSnap::RemoveSelfLoops(Graph);
    NumNodes = Graph->GetNodes();
    NumEdges = Graph->GetEdges();
    StatV = TFltV(fsLast);
    MotV = TFltV(fmLast);
  }
  void SetFracDeg();                                                 //Computes Fraction of Nodes with Degree >= AvgDeg and >= 2 * AvgDeg
  void SetBiCon();                                                   //Computes Fraction of Nodes in Largest Biconnected Component
  void SetDegCentr();                                                //Computes Degree Centrality values
  void SetBetwCentr();                                               //Computes Betweeness Centrality values
  void SetClsCentr();                                                //Computes Closeness Centrality values
  void SetKCore();                                                   //Computes k-core values
  void SetCCF();                                                     //Computes global and local clustering values
  void SetDiam();                                                    //Computes estimates of Diameter and Effective Diameter
  void SetFracLCC();                                                 //Computes fraction of nodes in largest and 2nd largest connected component
  void SetEgnVal();                                                  //Computes top 2 eigenvalues
  void SetAssty();                                                   //Computes assortativity coefficient
  void SetMotif();                                                   //Compute Motif Values
  void SetAll();                                                     //For computing everything at one go
};

