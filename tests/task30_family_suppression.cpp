#include "../src/merge_r8qkeep/family.cpp"
#include <iostream>
#include <stdexcept>

// Exercise actual family enumeration, including prior-part suppression.
// Index-only tests cannot detect a pair suppressed by a part that missed it.
int main(){
    for(const char* radius:{"1","2"}){
        if(setenv("FQDUP_FAMILY_RADIUS",radius,1))throw std::runtime_error("setenv failed");
        for(int threads:{1,8})for(int second_part:{0,2})for(bool reverse:{false,true}){
            std::string a(40,'A'),b=a;
            b[0]=b[1]='C';
            if(second_part)b[13]=b[14]='C';
            // N mismatch handling is asymmetric: repeat in both query orders.
            for(int npos:{-1,0,13,26}){
                auto left=a,right=b;
                if(npos>=0)left[npos]='N';
                std::vector<std::string> keys=reverse?std::vector<std::string>{right,left}:std::vector<std::string>{left,right};
                std::vector<int> parent{0,1};
                if(!neighbor_families(keys,parent,threads))throw std::runtime_error("unsupported test key");
                auto root=[&](int i){while(parent[i]!=i)i=parent[i];return i;};
                if(root(0)!=root(1))throw std::runtime_error("radius suppression lost a <=5 mismatch pair");
            }
        }
    }
    std::cout<<"TASK30_SUPPRESSION_PASS\n";
}
