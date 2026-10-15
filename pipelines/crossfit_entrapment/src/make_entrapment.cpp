// Paired 1:1 entrapment database for an already FIXED search database (A-16 / D-R03).
// Derived from pipelines/cohort_full_route_controls/make_controls.cpp (A91): the functions fnv, next_rng,
// cuts, shuf, validate and hexid are copied unchanged. Changes: no 1% hash sampling (every eligible protein
// gets one control), all original records are written first (sequence on one line, header verbatim), then
// the controls; the map records control -> source accession.
// Control construction (unchanged from A91): residues permuted within every zero-missed-cleavage tryptic
// block, K/R/P positions kept, per-protein seed fnv(sequence) ^ seed, so identical source sequences in two
// databases (e.g. FL and MG) receive identical controls. Controls equal to the source after I->L, proteins
// shorter than 9 residues and proteins with non-standard residues get no control (counted).
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
 if(argc!=6)throw runtime_error("usage: make_entrapment target.fasta combined.fasta map.tsv stats.json seed");
 ifstream in(argv[1]);ofstream out(argv[2]),map(argv[3]);if(!in||!out||!map)throw runtime_error("input/output unavailable");
 uint64_t seed=stoull(argv[5]),total=0,shortp=0,written=0,unchanged=0,duplicate=0,ambiguous=0;
 unordered_set<string> seen;vector<pair<string,string>> ctl;
 map<<"accession\tsource_accession\tlength\tconstruction_seed\n";
 string a,h,s,line;
 auto emit=[&](){if(h.empty())return;total++;out<<h<<"\n"<<s<<"\n";
   if(s.find_first_not_of("ACDEFGHIKLMNPQRSTVWY")!=string::npos){ambiguous++;return;}if(s.size()<9){shortp++;return;}
   uint64_t hh=fnv(s);
   string q;for(uint64_t attempt=0;attempt<4;attempt++){q=shuf(s,hh^seed^(attempt*0x9e3779b97f4a7c15ULL));string x=s,y=q;replace(x.begin(),x.end(),'I','L');replace(y.begin(),y.end(),'I','L');if(x!=y)break;}
   string x=s,y=q;replace(x.begin(),x.end(),'I','L');replace(y.begin(),y.end(),'I','L');if(x==y){unchanged++;return;}
   validate(s,q);string id="FLX"+hexid(q);if(!seen.insert(id).second){duplicate++;return;}
   ctl.push_back({id,q});map<<id<<"\t"<<a<<"\t"<<q.size()<<"\t"<<seed<<"\n";written++;
 };
 while(getline(in,line)){if(!line.empty()&&line.back()=='\r')line.pop_back();if(line.empty())continue;if(line[0]=='>'){emit();h=line;a=line.substr(1);a=a.substr(0,a.find_first_of(" \t"));if(a.rfind("FLX",0)==0)throw runtime_error("target accession already uses the FLX control prefix");s.clear();}else s+=line;}emit();
 for(auto&c:ctl)out<<">"<<c.first<<"\n"<<c.second<<"\n";
 if(in.bad()||!out||!map)throw runtime_error("stream failure");out.close();map.close();ofstream stats(argv[4]);
 stats<<"{\"target_records\":"<<total<<",\"controls_written\":"<<written<<",\"ambiguous_no_control\":"<<ambiguous<<",\"short_no_control\":"<<shortp<<",\"unchanged_after_il_no_control\":"<<unchanged<<",\"duplicate_control_ids_skipped\":"<<duplicate<<",\"seed\":"<<seed<<",\"all_written_cleavage_coordinates_verified\":true,\"all_written_block_compositions_verified\":true}\n";
 cerr<<"targets "<<total<<"; controls "<<written<<"\n";return written?0:2;
 }catch(const exception&e){cerr<<e.what()<<"\n";return 1;}}
