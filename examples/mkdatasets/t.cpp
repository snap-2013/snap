#include <iostream>
#include <cstdio>
#include <cstring>
#include <vector>
#include <algorithm>
#define MAXN 6000 
using namespace std;

int N, tri, path;
bool AM[MAXN][MAXN]; /// DAM is for directed graph
vector <int> G[MAXN], V; //D is the directed graph

int main(){
  int a,b;
  while(scanf("%d%d",&a,&b) == 2){
    AM[a][b] = AM[b][a] = 1;
    G[a].push_back(b);
    G[b].push_back(a);    
    N = max(N, max(a,b));
  }
  printf("N = %d\n",N);
  for(int i=0;i<=N;++i){
    if(G[i].empty()) continue;
    //Remove duplicates
    sort(G[i].begin(),G[i].end());
    V.clear();
    for(int j=0;j<(int)G[i].size();++j)
      if(j == 0 || G[i][j] != G[i][j-1]) V.push_back(G[i][j]);
    G[i].clear();
    for(int j=0;j<(int)V.size();++j) G[i].push_back(V[j]);

    for(int j=0;j<(int)G[i].size();++j)
      for(int k=j+1;k<(int)G[i].size();++k){
        int b = G[i][j], c = G[i][k];
        if(i != b && i != c)
          if(AM[b][c]) ++tri;
          else ++path;
      }
  }
  

  //Compute Clustering Coefficient
  double sum_coe = 0.0;

  for(int i=0;i<=N;++i){
    if(G[i].size() <= 1) continue;
    int tmp = 0;  
    for(vector<int>::iterator a=G[i].begin();a!=G[i].end();++a)
      for(vector<int>::iterator b=G[i].begin();b!=G[i].end();++b)
        if(AM[*a][*b]) ++tmp;
    //printf("tmp = %d, size = %d\n",G[i].size() * (G[i].size() - 1));
    //cin.get();
    double ret = (double)(tmp) / (double)(G[i].size() * (G[i].size() - 1)); 
   sum_coe += ret;
  }
 printf("triangles = (%d, %d), paths = %d\n",tri,tri/3,path);
  printf("Average clusterin coefficient is %.5lf / %d\n",sum_coe, N);
}
