// Streaming, composition/cleavage-preserving linked protein controls.
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;
uint64_t next64(uint64_t &x){x+=0x9e3779b97f4a7c15ULL;uint64_t z=x;z=(z^(z>>30))*0xbf58476d1ce4e5b9ULL;z=(z^(z>>27))*0x94d049bb133111ebULL;return z^(z>>31);}
string control(string s,uint64_t seed){
 for(char &c:s)if(c=='I')c='L';
 string out=s;size_t start=0;
 for(size_t end=1;end<=s.size();++end){
  if(end<s.size() && !((s[end-1]=='K'||s[end-1]=='R') && s[end]!='P'))continue;
  vector<size_t> pos;uint64_t state=seed ^ (end-start);
  // Position-only permutation: no dependence on peptide residue order/content.
  // Preserve K/R/P positions and the first residue of each cleavage segment.
  for(size_t j=start;j<end;++j){
   bool fixed=j==start||s[j]=='K'||s[j]=='R'||s[j]=='P';
   state^=uint64_t(fixed?static_cast<unsigned char>(s[j]):0)+0x9e3779b97f4a7c15ULL+(state<<6)+(state>>2);
   if(!fixed)pos.push_back(j);
  }
  auto perm=pos;
  for(size_t j=perm.size();j>1;--j)swap(perm[j-1],perm[next64(state)%j]);
  for(size_t j=0;j<pos.size();++j)out[pos[j]]=s[perm[j]];
  start=end;
 }
 return out;
}
int main(int argc,char **argv){
 try{
  if(argc!=5)throw runtime_error("usage: build_fused INPUT OUTPUT SEED EXPECTED_RECORDS");
  ifstream in(argv[1]);ofstream out(argv[2]);if(!in||!out)throw runtime_error("open failed");
  uint64_t seed=stoull(argv[3]),expected=stoull(argv[4]),n=0,unchanged=0;string line,head,seq;
  auto emit=[&](){if(head.empty())return;string e=control(seq,seed),norm=seq;replace(norm.begin(),norm.end(),'I','L');unchanged+=e==norm;out<<head<<'\n'<<seq<<string(51,'X')<<e<<'\n';++n;if(n%10000000==0)cerr<<n<<" proteins\n";};
  while(getline(in,line)){if(!line.empty()&&line.back()=='\r')line.pop_back();if(!line.empty()&&line[0]=='>'){emit();head=line;seq.clear();}else seq+=line;}emit();out.flush();if(!out)throw runtime_error("write failed");if(n!=expected)throw runtime_error("record count mismatch");
  cout<<"{\"proteins\":"<<n<<",\"unchanged_controls\":"<<unchanged<<",\"seed\":"<<seed<<",\"separator_length\":51}\n";
 }catch(const exception &e){cerr<<e.what()<<'\n';return 1;}
}
