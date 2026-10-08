#include "../src/merge_r8qkeep/repair.cpp"
#include <random>
using namespace fqdup::r8qkeep;
static void require(bool b){if(!b)throw std::runtime_error("repair skip unit failure");}
static std::string prior(int component){std::ostringstream o;o<<"FQDUP_PAIR_MIXTURE_V1 64\n";for(int f=0;f<2;++f){o<<f<<' '<<(component>=0);for(int z=0;z<66;++z)o<<' '<<(z==component?1:0);o<<'\n';}return o.str();}
int main(){
    std::mt19937 rng(23);uint64_t skipped=0,changed=0;
    for(int component:{-1,0,1,32,65}){
        const auto p=prior(component);RepairPipeline pipeline(p,p);std::string input;
        for(int row=0;row<32;++row){
            std::string a(40,'A'),q(40,'I');for(char& c:a)c="ACGT"[rng()%4];
            auto b=fqdup_r8qkeep_repair_internal::fqalign::reverse_complement(a);
            if(row%3==0)b[19]=b[19]=='A'?'C':'A';
            auto seq=a,qual=q;if(row%5==0)qual[0]='!';if(row%7==0)seq[20]='N';
            if(row%11==0){a[3]='N';b[5]='N';}
            input+="@x\t"+a+"\t"+q+"\t@y\t"+b+"\t"+q+"\tmerged\t@x\t"+seq+"\t"+qual+"\t\t\t\n";
        }
        input+="@x\tA\tI\t@y\tT\tI\tdropped\t\t\t\t\t\t\n";
        input+="@x\tACGT\tIIII\t@y\tACGT\tIIII\torphan_r1\t@x\tACGT\tIIII\t\t\t\n";
        for(double cutoff:{0.0,.0001,.01,.316,1.0}){
            RepairShadow s;const auto reference=pipeline.repair(input,9,cutoff,0,1,&s);
            require(reference==pipeline.repair(input,9,cutoff));require(s.skip_but_changed==0);
            require(s.skip_and_unchanged+s.noskip_unchanged+s.noskip_changed==34);
            RepairShadow parallel;require(reference==pipeline.repair(input,9,cutoff,0,4,&parallel));
            require(s.skip_and_unchanged==parallel.skip_and_unchanged.load());
            require(s.noskip_changed==parallel.noskip_changed.load());
            std::vector<bool> excluded(34,true);RepairShadow excluded_counts;
            require(reference==pipeline.repair(input,9,cutoff,0,1,&excluded_counts,&excluded));
            require(excluded_counts.skip_and_unchanged+excluded_counts.noskip_unchanged+excluded_counts.noskip_changed==0);
            skipped+=s.skip_and_unchanged;changed+=s.noskip_changed;
        }
    }
    // Exact budget boundary and monotonicity when the nominal model is unsupported.
    const auto p=prior(-1);auto model=fqdup_r8qkeep_repair_internal::model(p);RepairSkipBound bound(model);
    auto v=fqdup_r8qkeep_repair_internal::split("@x\tA\tI\t@y\tT\tI\tmerged\t@x\tA\tI\t\t\t");
    require(!bound.unchanged(v,model,.000099999999));require(bound.unchanged(v,model,.0001));require(bound.unchanged(v,model,1));
    v[9]=" ";require(!bound.unchanged(v,model,1));v[9]="I";v[8]="AA";require(!bound.unchanged(v,model,1));
    // Exercise multiple in-flight chunks, nonzero global indices, and mixed exclusions.
    {
        RepairPipeline pipeline(p,p);std::string input;
        for(int i=0;i<1100;++i)input+="@x\tA\tI\t@y\tT\tI\tmerged\t@x\tA\tI\t\t\t\n";
        std::vector<bool> excluded(1107);for(size_t i=7;i<excluded.size();++i)excluded[i]=(i%2==0);
        RepairShadow serial,parallel;
        require(pipeline.repair(input,9,.316,7,1,&serial,&excluded)==pipeline.repair(input,9,.316,7,8,&parallel,&excluded));
        require(serial.skip_and_unchanged==550&&parallel.skip_and_unchanged==550);
        require(parallel.skip_but_changed==0&&parallel.noskip_changed==0);
        RepairShadow mode7;pipeline.repair(input,7,.316,7,8,&mode7);
        require(mode7.skip_and_unchanged+mode7.skip_but_changed+mode7.noskip_unchanged+mode7.noskip_changed==0);
    }
    // Stress the normalization/dot-product allowance independently of alignment.
    for(int trial=0;trial<10000;++trial){
        std::array<double,66> post{},cp{};double mass=0,maximum=0;
        for(int z=0;z<66;++z){post[z]=z==trial%66?1:std::exp(-double(rng()%746));mass+=post[z];
            const auto risk=fqdup_r8qkeep_repair_internal::fqquality::PairMixture::conditional(
                fqdup_r8qkeep_repair_internal::fqalign::phred_error(char(33+rng()%94)),
                fqdup_r8qkeep_repair_internal::fqalign::phred_error(char(33+rng()%94)),trial%2);
            cp[z]=trial%3?risk.first:risk.second;maximum=std::max(maximum,cp[z]);}
        double value=0;for(int z=0;z<66;++z){post[z]/=mass;value+=post[z]*cp[z];}
        require(value<=maximum+1e-12);
    }
    require(skipped>0&&changed>0);std::cout<<"TASK23_REPAIR_SKIP_UNIT_PASS skipped="<<skipped<<" changed="<<changed<<'\n';
}
