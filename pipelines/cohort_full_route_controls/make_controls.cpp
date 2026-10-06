// Reservoir-sampled synthetic controls, independent of cohort predictions.
// Preserve composition within every zero-missed-cleavage tryptic block, and
// every K/R/P position. Verify cleavage coordinates and block compositions.
// This is NOT a globally bijective peptide-pair construction or FDP estimator.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>
using namespace std;
uint64_t fnv(const string&s,uint64_t h=14695981039346656037ULL){for(unsigned char c:s){h^=c;h*=1099511628211ULL;}return h;}
uint64_t next_rng(uint64_t &x){x+=0x9e3779b97f4a7c15ULL;uint64_t z=x;z=(z^(z>>30))*0xbf58476d1ce4e5b9ULL;z=(z^(z>>27))*0x94d049bb133111ebULL;return z^(z>>31);}
vector<size_t> cuts(const string&s){vector<size_t>c{0};for(size_t i=0;i<s.size();i++)if((s[i]=='K'||s[i]=='R')&&(i+1==s.size()||s[i+1]!='P'))c.push_back(i+1);if(c.back()!=s.size())c.push_back(s.size());return c;}
string shuf(const string&s,uint64_t seed){string o=s;auto c=cuts(s);for(size_t k=1;k<c.size();k++){vector<size_t>p;for(size_t i=c[k-1];i<c[k];i++)if(s[i]!='K'&&s[i]!='R'&&s[i]!='P')p.push_back(i);for(size_t n=p.size();n>1;n--){size_t j=next_rng(seed)%n;swap(o[p[n-1]],o[p[j]]);}}return o;}
void validate(const string&a,const string&b){auto c=cuts(a);if(c!=cuts(b))throw runtime_error("cleavage changed");for(size_t k=1;k<c.size();k++){array<int,26>x{};for(size_t i=c[k-1];i<c[k];i++){x[a[i]-'A']++;x[b[i]-'A']--;}for(int n:x)if(n)throw runtime_error("block composition changed");}}
string hexid(const string&s){ostringstream o;o<<hex<<setfill('0')<<setw(16)<<fnv(s)<<setw(16)<<fnv(s,1099511628211ULL);return o.str();}
int main(int argc,char**argv){try{
 if(argc!=6)throw runtime_error("usage: make_controls reservoir.fasta controls.fasta map.tsv stats.json seed");
 ifstream in(argv[1]);ofstream out(argv[2]),map(argv[3]);if(!in||!out||!map)throw runtime_error("input/output unavailable");
 uint64_t seed=stoull(argv[5]),total=0,eligible=0,sampled=0,written=0,unchanged=0,duplicate=0,ambiguous=0;
 unordered_set<string> seen;
 map<<"accession\tclass\tsource_accession\tlength\tconstruction_seed\n";
 string a,s,line;
 auto emit=[&](){if(a.empty())return;total++;if(s.find_first_not_of("ACDEFGHIKLMNPQRSTVWY")!=string::npos){ambiguous++;return;}if(s.size()<9)return;eligible++;
   uint64_t h=fnv(s);if((h^seed)%100!=0)return;sampled++;
   string q;for(uint64_t attempt=0;attempt<4;attempt++){q=shuf(s,h^seed^(attempt*0x9e3779b97f4a7c15ULL));string x=s,y=q;replace(x.begin(),x.end(),'I','L');replace(y.begin(),y.end(),'I','L');if(x!=y)break;}
   string x=s,y=q;replace(x.begin(),x.end(),'I','L');replace(y.begin(),y.end(),'I','L');if(x==y){unchanged++;return;}
   validate(s,q);string id="FLX"+hexid(q);if(!seen.insert(id).second){duplicate++;return;}
   out<<">"<<id<<"\n"<<q<<"\n";map<<id<<"\tshuffle\t"<<a<<"\t"<<q.size()<<"\t"<<seed<<"\n";written++;
 };
 while(getline(in,line)){if(!line.empty()&&line.back()=='\r')line.pop_back();if(line.empty())continue;if(line[0]=='>'){emit();a=line.substr(1);a=a.substr(0,a.find_first_of(" \t"));s.clear();}else s+=line;}emit();
 if(in.bad()||!out||!map)throw runtime_error("stream failure");out.close();map.close();ofstream stats(argv[4]);
 stats<<"{\"input_sequences\":"<<total<<",\"eligible_standard_sequences\":"<<eligible<<",\"ambiguous_sequences_excluded\":"<<ambiguous<<",\"sampled_at_hash_fraction_0_01\":"<<sampled<<",\"written\":"<<written<<",\"unchanged_after_il_excluded\":"<<unchanged<<",\"duplicate_control_ids_excluded\":"<<duplicate<<",\"seed\":"<<seed<<",\"all_written_cleavage_coordinates_verified\":true,\"all_written_block_compositions_verified\":true,\"globally_bijective_peptide_pairs_claimed\":false}\n";
 cerr<<"Scanned "<<total<<"; wrote "<<written<<" verified control proteins\n";return written?0:2;
 }catch(const exception&e){cerr<<e.what()<<"\n";return 1;}}
