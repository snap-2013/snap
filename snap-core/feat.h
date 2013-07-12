#define THRESHOLD_VAL 30000

typedef enum { Mn=0, Md, Mx, Avg, Std, EnumStatLast } EnumStat;
typedef enum { First=0, Second, EnumSizeLast } EnumSize;

class TGraphFeature {
private:
  PUNGraph Graph;
public:
  TInt NumNodes, NumEdges;
  int64 OpenTriads, ClosedTriads; TFlt FracTriads;       //Triads and Triangles
  TFlt CCF_Global; 
  TTuple <TFlt, EnumStatLast> CCF;                                  //Clustering
  TTuple <TFlt, EnumStatLast> DegC;                                 //Degree Centrality
  TTuple <TFlt, EnumStatLast> BetwC;                                //Betweeness Centrality
  TTuple <TFlt, EnumStatLast> ClsC;                                 //Closeness Centrality
  TInt FullDiam;   TFlt EffDiam;                         //Distance
  TFlt A, B, SigA, SigB, Chi2, R2;                       //PowerFit Variables
  TTuple <TFlt, EnumSizeLast> FracLCC;                              //Fraction of vertices in Largest 2 Connected Components
  TInt KCoreSize; TInt FracKCore;                        //K-Core values
  TTuple <TFlt, EnumSizeLast> EgnVal;                               //Largest 2 Eigenvalues
  TFlt AsstyCor;                                         //Assortativity Coefficient
  TTuple <TFlt, 10> motif;                               //Array containing fraction of each motif occuring
private:
  void GrabStats (TFltV& V, TTuple <TFlt, EnumStatLast> & Stats);
public:
  TGraphFeature (PUNGraph& _Graph) : Graph(_Graph) {
    TSnap::RenumberNodes(Graph);
    TSnap::RemoveSelfLoops(Graph);
    NumNodes = Graph->GetNodes();
    NumEdges = Graph->GetEdges();
  }
  void SetDegPowerFit();
  void SetDegCentr();
  void SetBetwCentr();
  void SetClsCentr();
  void SetKCore();
  void SetCCF();
  void SetDiam();
  void SetFracLCC();
  void SetEgnVal();
  void SetAssty();
  void SetAll();                                         //For computing everything at one go
};

