#include "../src/merge_r8qkeep/repair.cpp"
#include <random>
using namespace fqdup::r8qkeep;
static void require(bool ok){if(!ok)throw std::runtime_error("clean differential failure");}
static std::string prior(int component){
    std::ostringstream out;out<<"FQDUP_PAIR_MIXTURE_V1 64\n";
    for(int f=0;f<2;++f){out<<f<<' '<<(component>=0);for(int j=0;j<66;++j)out<<' '<<(j==component?1:0);out<<'\n';}
    return out.str();
}
int main(){
    // Cached conditional sums must match the original order bit-for-bit,
    // including collisions, both agreement states and subnormal posteriors.
    for(int variant=0;variant<3;++variant){
        std::array<double,fqquality::PairMixture::K> post{};
        for(int j=0;j<int(post.size());++j)post[j]=variant==0?1.0/post.size():variant==1?std::ldexp(1.0,-j-1000):(j==0?1.0:0.0);
        fqquality::PairMixture::ConditionalProducts products(post);
        for(int repeat=0;repeat<2;++repeat)for(int q1=33;q1<=126;q1+=7)for(int q2=33;q2<=126;q2+=11)for(bool agree:{false,true}){
            double p1=0,p2=0;
            for(int z=0;z<int(post.size());++z){auto risk=fqquality::PairMixture::lookup(z,char(q1),char(q2),agree).risk;p1+=post[z]*risk.first;p2+=post[z]*risk.second;}
            const auto expected=std::make_pair(p1,p2);
            require(products.get(char(q1),char(q2),agree)==expected);
            require(products.get(char(q1),char(q2),agree)==expected);
        }
    }

    std::mt19937 rng(2503);std::size_t positives=0,negatives=0;
    for(int component:{-1,0,40,65}){
        RepairPipeline pipeline(prior(component),prior(component));
        for(int trial=0;trial<60;++trial){
            const int insert=35+trial,n=trial%2?insert-4:insert+7,m=trial%2?insert-5:insert+9;
            std::string molecule(insert,'A');for(char& c:molecule)c="ACGT"[rng()%4];
            std::string a=molecule.substr(0,n),b=fqalign::reverse_complement(molecule.substr(std::max(0,insert-m)));
            while(int(a.size())<n)a+="ACGT"[rng()%4];
            while(int(b.size())<m)b+="ACGT"[rng()%4];
            std::string q1(n,char(33+(trial%3?40:8))),q2(m,char(33+(trial%3?40:8)));
            std::vector<std::string> v{"@pair/1",a,q1,"@pair/2",b,q2,"merged","@pair/1",molecule,std::string(insert,'I'),"","",""};
            if(trial%6==1)v[4][15]=v[4][15]=='A'?'C':'A';
            if(trial%6==2)v[1][15]='N';
            if(trial%6==3)v[8][15]='N';
            if(trial%6==4){v[6]="dropped";for(int j=7;j<13;++j)v[j].clear();}
            if(trial%6==5){v[6]="unmerged";v[7]=v[0];v[8]=v[1];v[9]=v[2];v[10]=v[3];v[11]=v[4];v[12]=v[5];}
            std::string text;for(std::size_t j=0;j<v.size();++j){if(j)text+='\t';text+=v[j];}text+='\n';
            RepairAlignmentCache cache(0,1);
            const auto reference=pipeline.route(text,true),fast=pipeline.route(text,false,1,0,&cache);
            require(reference.routed==fast.routed&&reference.replacements==fast.replacements);
            for(double cut:{0.0,.0001,.001,.01,.1,.316,.75,1.0}){
                RepairShadow shadow;
                const auto expected=pipeline.repair(text,9,cut,0,1,&shadow,&reference.routed,&reference.repair_bounds,nullptr,&fast.clean_bounds);
                const auto actual=pipeline.repair(text,9,cut,0,1,nullptr,&fast.routed,&fast.repair_bounds,&cache,&fast.clean_bounds);
                require(expected==actual&&shadow.clean_but_changed==0&&shadow.skip_but_changed==0);
                positives+=shadow.clean_and_unchanged;negatives+=shadow.notclean;
                require(shadow.clean_and_unchanged+shadow.notclean==1);
                const auto selected=pipeline.repair(text,7,cut,0,1,nullptr,nullptr,nullptr,&cache,nullptr,&fast.routed);
                if(fast.routed[0])require(selected==pipeline.repair(text,7,cut));
                else require(selected=="0\tdropped\t\t\t\t\t\t\n");
            }
        }
    }
    require(positives>0&&negatives>0);
    std::cout<<"TASK25_CLEAN_PASS clean="<<positives<<" notclean="<<negatives<<'\n';
}
