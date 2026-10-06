#include <fstream>
#include <iostream>
#include <vector>
#include <array>
#include <string>
#include <stdexcept>
#include <cstdint>
int main(int argc,char**argv){
 try{
  if(argc!=3)throw std::runtime_error("usage: composition input.fasta output_prefix");
  std::ifstream f(argv[1],std::ios::binary);
  if(!f)throw std::runtime_error("cannot read FASTA");
  std::array<uint64_t,26> aa{}; std::vector<uint64_t> hist(1024,0);
  uint64_t entries=0,residues=0,length=0,stops=0;bool started=false;
  auto finish=[&](){
    if(!started)return;
    if(length==0)throw std::runtime_error("empty FASTA entry");
    if(length>=hist.size())hist.resize(length+1,0);
    hist[length]++;entries++;residues+=length;
  };
  std::string line;
  while(std::getline(f,line)){
   if(!line.empty()&&line.back()=='\r')line.pop_back();
   if(line.empty())continue;
   if(line[0]=='>'){finish();started=true;length=0;}
   else{
    if(!started)throw std::runtime_error("sequence before header");
    for(unsigned char c:line){
     if(c>='A'&&c<='Z')aa[c-'A']++;
     else if(c=='*')stops++;
     else throw std::runtime_error("unexpected residue byte");
     length++;
    }
   }
  }
  if(f.bad())throw std::runtime_error("FASTA read failure");
  finish();
  std::ofstream o(std::string(argv[2])+".json");
  o<<"{\n\"entries\":"<<entries<<",\n\"residues\":"<<residues<<",\n\"stop_symbols\":"<<stops<<",\n\"aa_counts\":{";
  for(int i=0;i<26;i++){if(i)o<<",";o<<"\""<<char('A'+i)<<"\":"<<aa[i];}
  o<<"},\n\"length_histogram\":{";bool comma=false;
  for(size_t i=0;i<hist.size();i++)if(hist[i]){if(comma)o<<",";comma=true;o<<"\""<<i<<"\":"<<hist[i];}
  o<<"}\n}\n";o.close();if(!o)throw std::runtime_error("output write failure");
  std::cout<<"entries="<<entries<<" residues="<<residues<<"\n";
 }catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}
}

