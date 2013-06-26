#include <iostream>
#include <cstdio>
#include <map>
#include <vector>
using namespace std;

vector <pair<int,int> > E;
map<int,int> M;

int main(){
  int a,b;
  while(scanf("%d%d",&a,&b) == 2) M[a] = M[b] = 1, E.push_back(pair<int,int>(a,b));
  int cnt = 0;
  for(map<int,int>::iterator it=M.begin();it!=M.end();++it) it->second = cnt++;
  printf("cnt = %d\n",cnt);
  freopen("b.in","w",stdout);
  for(vector<pair<int,int> >::iterator it=E.begin();it!=E.end();++it)
    printf("%d %d\n",M[it->first],M[it->second]);
}
