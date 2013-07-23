namespace TSnap {

/*
 * For MotifSize = 3, 0 is ClosedTriads, 1 is OpenTriads
 *
 * For MotifSize = 4, we use the following convention to label the motifs
 *
 * NOTE: The pair (a,b) represents an edge between a and b 
 *
 * 0 : (1,2), (2,3), (3,4)
 *
 * 1 : (1,2), (2,3), (3,4), (1,4)
 *
 * 2 : (1,2), (1,3), (1,4)
 *
 * 3 : (1,2), (1,3), (1,4), (2,3)
 *
 * 4 : (1,2), (1,3), (1,4), (2,3), (3,4)
 *
 * 5 : (1,2), (1,3), (1,4), (2,3), (3,4), (2,4)
 *
 */

void GetMotifCount(const PUNGraph& G, const int MotifSize, TVec <int64> & MotifV, const int num);

}; //namespace TSnap
