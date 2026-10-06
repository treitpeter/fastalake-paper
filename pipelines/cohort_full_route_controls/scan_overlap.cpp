// Exact I/L-normalised query presence in the ORIGINAL reservoir, across all
// protein substrings. Protein boundaries reset the automaton; line wraps do not.
#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <queue>
#include <stdexcept>
#include <string>
#include <vector>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <omp.h>
using namespace std;
struct Node{array<int,26>go{};int fail=0;vector<int>out;};
int main(int argc,char**argv){try{
 if(argc!=5)throw runtime_error("usage: scan_overlap queries.txt reservoir.fasta output.tsv threads");
 vector<Node>trie(1);vector<string>queries;string p;ifstream in(argv[1]);if(!in)throw runtime_error("queries unavailable");
 while(getline(in,p)){if(p.empty())continue;replace(p.begin(),p.end(),'I','L');int state=0;for(char c:p){if(c<'A'||c>'Z')throw runtime_error("bad query");int n=trie[state].go[c-'A'];if(!n){n=trie.size();trie[state].go[c-'A']=n;trie.emplace_back();}state=n;}trie[state].out.push_back(queries.size());queries.push_back(p);}
 queue<int>q;for(int c=0;c<26;c++)if(trie[0].go[c])q.push(trie[0].go[c]);
 while(!q.empty()){int v=q.front();q.pop();int f=trie[v].fail;trie[v].out.insert(trie[v].out.end(),trie[f].out.begin(),trie[f].out.end());for(int c=0;c<26;c++){int w=trie[v].go[c];if(w){trie[w].fail=trie[f].go[c];q.push(w);}else trie[v].go[c]=trie[f].go[c];}}
 int fd=open(argv[2],O_RDONLY);if(fd<0)throw runtime_error("reservoir unavailable");struct stat st;if(fstat(fd,&st))throw runtime_error("stat failed");size_t n=st.st_size;char*data=(char*)mmap(nullptr,n,PROT_READ,MAP_PRIVATE,fd,0);if(data==MAP_FAILED)throw runtime_error("mmap failed");
 int threads=stoi(argv[4]);vector<vector<unsigned char>>seen(threads,vector<unsigned char>(queries.size(),0));
 #pragma omp parallel num_threads(threads)
 {
  int t=omp_get_thread_num(), nt=omp_get_num_threads();size_t start=n*t/nt,end=n*(t+1)/nt;
  if(start)while(start<n&&!(data[start]=='>'&&(start==0||data[start-1]=='\n')))start++;
  if(end<n)while(end<n&&!(data[end]=='>'&&(end==0||data[end-1]=='\n')))end++;
  int state=0;bool header=false;
  for(size_t i=start;i<end;i++){char c=data[i];if(c=='>'){header=true;state=0;continue;}if(header){if(c=='\n')header=false;continue;}if(c=='\n'||c=='\r')continue;if(c=='I')c='L';if(c<'A'||c>'Z'){state=0;continue;}state=trie[state].go[c-'A'];for(int id:trie[state].out)seen[t][id]=1;}
 }
 ofstream out(argv[3]);out<<"peptide\tpresent_in_original_reservoir\n";size_t found=0;for(size_t i=0;i<queries.size();i++){bool hit=false;for(int t=0;t<threads;t++)hit|=seen[t][i];out<<queries[i]<<'\t'<<int(hit)<<'\n';found+=hit;}munmap(data,n);close(fd);cerr<<"Queries "<<queries.size()<<"; present "<<found<<"; states "<<trie.size()<<"\n";
 }catch(const exception&e){cerr<<e.what()<<'\n';return 1;}return 0;}
