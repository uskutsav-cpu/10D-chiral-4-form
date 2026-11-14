// Exact Sym^n(126) character, compressed only by Weyl D5 weight orbits.
// No graph ranks or target Hilbert coefficients are inputs.
#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>
using Weight=std::array<int,5>;
using Integer=__int128_t;
using Character=std::unordered_map<uint32_t,Integer>;
std::string decimal(Integer n){if(n==0)return "0";bool neg=n<0;if(neg)n=-n;std::string s;while(n){s.push_back('0'+n%10);n/=10;}if(neg)s.push_back('-');std::reverse(s.begin(),s.end());return s;}
uint32_t key(const Weight& a){uint32_t k=0;for(int i=0;i<4;i++)k=(k<<6)|a[i];return (k<<6)|(a[4]+32);}
Weight dominant(Weight a){int sign=1;for(int& x:a){if(x<0){sign=-sign;x=-x;}}std::sort(a.begin(),a.end(),std::greater<int>());a[4]*=sign;return a;}
int factorial(int n){int r=1;for(int i=2;i<=n;i++)r*=i;return r;}
int orbit_size(Weight a){for(int& x:a)x=std::abs(x);int z=std::count(a.begin(),a.end(),0);std::map<int,int> counts;for(int x:a)counts[x]++;int denominator=1;for(auto [x,n]:counts)denominator*=factorial(n);return factorial(5)*(1<<(z?5-z:4))/denominator;}
Integer binomial(int n,int k){Integer r=1;for(int j=1;j<=k;j++){r*=n-k+j;assert(r%j==0);r/=j;}return r;}
std::vector<std::pair<Weight,int>> weights(){std::vector<std::pair<Weight,int>> out;for(int code=0;code<243;code++){int q=code,support=0,product=1;Weight w;for(int& x:w){x=q%3-1;q/=3;support+=(x!=0);if(x)product*=x;}if(support==1)out.push_back({w,3});if(support==3||(support==5&&product==1))out.push_back({w,1});}int dim=0;for(auto [w,m]:out)dim+=m;assert(dim==126&&out.size()==106);return out;}
std::map<uint32_t,int> denominator(){std::map<uint32_t,int> out;Weight permutation={0,1,2,3,4},rho={4,3,2,1,0};int terms=0;do{int parity=1;for(int i=0;i<5;i++)for(int j=i+1;j<5;j++)if(permutation[i]>permutation[j])parity=-parity;for(int bits=0;bits<32;bits++){if(__builtin_popcount((unsigned)bits)%2)continue;Weight shift;for(int i=0;i<5;i++)shift[i]=rho[i]-((bits>>i&1)?-1:1)*rho[permutation[i]];out[key(dominant(shift))]+=parity;terms++;}}while(std::next_permutation(permutation.begin(),permutation.end()));assert(terms==1920);return out;}
int main(int argc,char**argv){int maximum=argc>1?std::stoi(argv[1]):22;if(maximum<0||maximum>22)return 2;std::string prefix=argc>2?argv[2]:"character";auto w=weights();auto den=denominator();std::vector<Character> h(maximum+1);h[0][key({0,0,0,0,0})]=1;std::ofstream data(prefix+"-weights.tsv"),summary(prefix+".json");summary<<"{\"schema\":1,\"representation\":\"D5 highest weight (1,1,1,1,1), Dynkin (0,0,0,0,2)\",\"degrees\":[\n";for(int n=0;n<=maximum;n++){auto start=std::chrono::steady_clock::now();if(n){for(int a=0;a<=n;a++)for(int b=0;b<=a;b++)for(int c=0;c<=b;c++)for(int d=0;d<=c;d++)for(int e=-d;e<=d;e++){if(((a+b+c+d+e-n)%2)!=0)continue;Weight lambda={a,b,c,d,e};Integer sum=0;for(int k=1;k<=n;k++){int m=n-k;for(auto [nu,multiplicity]:w){Weight shifted;bool possible=true;for(int i=0;i<5;i++){shifted[i]=lambda[i]-k*nu[i];if(std::abs(shifted[i])>m){possible=false;break;}}if(!possible)continue;auto it=h[m].find(key(dominant(shifted)));if(it!=h[m].end())sum+=multiplicity*it->second;}}assert(sum>=0&&sum%n==0);if(sum)h[n][key(lambda)]=sum/n;}}
Integer dimension=0, invariant=0;std::vector<std::pair<uint32_t,Integer>> ordered(h[n].begin(),h[n].end());std::sort(ordered.begin(),ordered.end());for(auto [code,m]:ordered){Weight lambda;uint32_t q=code;lambda[4]=int(q&63)-32;q>>=6;for(int i=3;i>=0;i--){lambda[i]=q&63;q>>=6;}dimension+=m*orbit_size(lambda);data<<n;for(int x:lambda)data<<'\t'<<x;data<<'\t'<<decimal(m)<<'\n';}assert(dimension==binomial(125+n,n));for(auto [code,sign]:den){auto it=h[n].find(code);if(it!=h[n].end())invariant+=sign*it->second;}assert(invariant>=0);if(n%2)assert(invariant==0);double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();if(n)summary<<",\n";summary<<"{\"degree\":"<<n<<",\"invariants\":\""<<decimal(invariant)<<"\",\"symmetric_dimension\":\""<<decimal(dimension)<<"\",\"nonzero_weight_orbits\":"<<h[n].size()<<",\"seconds\":"<<seconds<<"}";summary.flush();data.flush();std::cout<<"degree "<<n<<" invariants "<<decimal(invariant)<<" orbits "<<h[n].size()<<" seconds "<<seconds<<std::endl;}summary<<"\n]}\n";}
