const int MaxSamplingSize = 30000;

typedef enum { fsMn=0, fsMd, fsMx, fsAvg, fsStd, fsEnumStatLast } EnumStat;
typedef enum { fsFirst=0, fsSecond, fsEnumSizeLast } EnumSize;

class TGraphFeature {
private:
  PUNGraph Graph;
public:
  TInt NumNodes, NumEdges;
  int64 OpenTriads, ClosedTriads; TFlt FracTriads;                   //Triads and Triangles
  TFlt CCF_Global; 
  TTuple <TFlt, fsEnumStatLast> CCF;                                 //Clustering
  TTuple <TFlt, fsEnumStatLast> DegC;                                //Degree Centrality
  TTuple <TFlt, fsEnumStatLast> BetwC;                               //Betweeness Centrality
  TTuple <TFlt, fsEnumStatLast> ClsC;                                //Closeness Centrality
  TInt FullDiam;   TFlt EffDiam;                                     //Distance
  TFlt A, B, SigA, SigB, Chi2, R2;                                   //PowerFit Variables
  TTuple <TFlt, fsEnumSizeLast> FracLCC;                             //Fraction of vertices in Largest 2 Connected Components
  TInt KCoreSize; TFlt FracKCore;                                    //K-Core values
  TTuple <TFlt, fsEnumSizeLast> EgnVal;                              //Largest 2 Eigenvalues
  TFlt AsstyCor;                                                     //Assortativity Coefficient
  TTuple <TFlt, 10> motif;                                           //Array containing fraction of each motif occuring
private:
  void GrabStats (TFltV& V, TTuple <TFlt, fsEnumStatLast> & Stats);  //Updates the fsMn, fsMd, fsMx, fsAvg, fsStd values for the centrality values 
public:
  TGraphFeature (PUNGraph& _Graph) : Graph(_Graph) {                 //Constructor. Nodes relabelled from 0 to N-1 and self-loops removed.
    TSnap::RenumberNodes(Graph);                                     
    TSnap::RemoveSelfLoops(Graph);
    NumNodes = Graph->GetNodes();
    NumEdges = Graph->GetEdges();
  }
  void SetDegPowerFit();                                             //Computes the powerfit variables
  void SetDegCentr();                                                //Computes Degree Centrality values
  void SetBetwCentr();                                               //Computes Betweeness Centrality values
  void SetClsCentr();                                                //Computes Closeness Centrality values
  void SetKCore();                                                   //Computes k-core values
  void SetCCF();                                                     //Computes global and local clustering values
  void SetDiam();                                                    //Computes estimates of Diameter and Effective Diameter
  void SetFracLCC();                                                 //Computes fraction of nodes in largest and 2nd largest connected component
  void SetEgnVal();                                                  //Computes top 2 eigenvalues
  void SetAssty();                                                   //Computes assortativity coefficient
  void SetAll();                                                     //For computing everything at one go
};

