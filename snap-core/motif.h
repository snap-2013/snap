typedef enum {
  mtUndef = -1,
  mtThreeClosed, mtThreeOpen,
  mtLast 
}ThreeMotif;

typedef enum {
  mfUndef = -1,
  mfFourLine, mfFourTriangleEdge, mfFourStar, mfFourSquare, mfFourSquareDiag, mfFourComplete,
  mfLast
}FourMotif;

namespace TSnap {

/// Get Motif Counts for n = 3 or n = 4
void GetMotifCount(const PUNGraph& G, const int MotifSize, TVec <int64> & MotifV, const int num);

}; //namespace TSnap
