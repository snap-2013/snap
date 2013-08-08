typedef enum {
  fsUndef=-1,
  fsMnCCF, fsMdCCF, fsMxCCF, fsAvgCCF, fsStdCCF,
  fsMnDegC, fsMdDegC, fsMxDegC, fsAvgDegC, fsStdDegC,
  fsMnBetwC, fsMdBetwC, fsMxBetwC, fsAvgBetwC, fsStdBetwC,
  fsMnClsC, fsMdClsC, fsMxClsC, fsAvgClsC, fsStdClsC,
  fsFullDiam, fsEffDiam,
  fsFracLCC1, fsFracLCC2, fsFracDeg, fsFrac2Deg, fsFracBiCon,
  fsKCoreSize, fsFracKCore, fsFrac2Core, fsFrac3Core, fsFrac4Core,
  fsEgnVal1, fsEgnVal2, fsAsstyCor,
  fsLast
}TGraphFeatureStats;

typedef enum {
  fmUndef=-1,
  fmThreeClosed, fmThreeOpen,
  fmFourLine, fmSquare, fmFourStar, fmFourTriangleEdge, fmFourSquareDiag, fmFourComplete,
  fmLast
}TGraphFeatureMotifs;

const int MaxSamplingSize = 10000;
const int MaxMotSampleSize = 20000;

const TStr fsName[fsLast] = {
  "fsMnCCF","fsMdCCF","fsMxCCF","fsAvgCCF","fsStdCCF",
  "fsMnDegC","fsMdDegC","fsMxDegC","fsAvgDegC","fsStdDegC",
  "fsMnBetwC","fsMdBetwC","fsMxBetwC","fsAvgBetwC","fsStdBetwC",
  "fsMnClsC","fsMdClsC","fsMxClsC","fsAvgClsC","fsStdClsC",
  "fsFullDiam","fsEffDiam",
  "fsFracLCC1","fsFracLCC2","fsFracDeg","fsFrac2Deg","fsFracBiCon",
  "fsKCoreSize","fsFracKCore","fsFrac2Core","fsFrac3Core","fsFrac4Core",
  "fsEgnVal1","fsEgnVal2","fsAsstyCor"
};

const TStr fmName[fmLast] = {
  "fmThreeClosed","fmThreeOpen",
  "fmFourLine","fmFourSquare","fmFourStar","fmFourTriangleEdge","fmFourSquareDiag","fmFourComplete"
};

class TGraphFeature {
  private:
    PUNGraph Graph;
  public:
    TInt NumNodes, NumEdges;
    TFltV StatV, MotV;
  private:
    /// Computes Minimum, Median, Maximum, Average and Standard Deviation for Local Centrality Values
    void SetStats (TFltV& V, TFlt & mn, TFlt & md, TFlt & mx, TFlt & avg, TFlt & std);    
  public:
    /// Constructor. Nodes relabelled from 0 to N-1 and Self-Loops removed
    TGraphFeature (PUNGraph& _Graph) : Graph(_Graph) {
      TSnap::RenumberNodes(Graph);                                     
      TSnap::RemoveSelfLoops(Graph);
      NumNodes = Graph->GetNodes();
      NumEdges = Graph->GetEdges();
      StatV = TFltV(fsLast);
      MotV = TFltV(fmLast);
    }
    /// Computes Fraction of Nodes with Degree >= AvgDeg and >= 2 * AvgDeg
    void SetFracDeg();                         
    /// Computes Fraction of Nodes in Largest Biconnected Component
    void SetBiCon();          
    /// Computes Degree Centrality Values
    void SetDegCentr();                   
    /// Computes Betweeness Centrality Values
    void SetBetwCentr();                   
    /// Computes Closeness Centrality Values
    void SetClsCentr();                   
    /// Computes Global and Local Clustering Values
    void SetCCF();    
    /// Computes Fraction of Nodes in Largest K-Core and in {2,3,4}-Core. 
    void SetKCore();                                                                                  
    /// Computes Estimates of Diameter and Effective Diameter
    void SetDiam();                                             
    /// Computes Fraction of Nodes in Largest and 2nd Largest Connected Component
    void SetFracLCC();                                        
    /// Computes Top 2 EigenValues
    void SetEgnVal();                               
    /// Computes Assortativity Coefficient
    void SetAssty();                             
    /// Computes Motif Values
    void SetMotif();              
    /// Computes all Values at one go
    void SetAll();
    /// Outputs Computed Values in Either Human-Friendly (0) or Machine-Friendly Format (1)
    void PrintVals(const bool format = 0); 
};

