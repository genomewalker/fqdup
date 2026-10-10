#include "../src/merge_r8qkeep/worker_placement.hpp"
#include <iostream>
#include <stdexcept>
#include <thread>
using namespace fqdup::r8qkeep;
static void require(bool ok){if(!ok)throw std::runtime_error("worker placement");}
int main(){
 const std::vector<WorkerCpu> two_sockets{{20,0,0},{22,0,1},{47,1,0},{116,0,0},{118,0,1},{143,1,0}};
 require(select_worker_cpus(two_sockets,3)==std::vector<int>({20,22,47}));
 require(select_worker_cpus(two_sockets,1).empty());require(select_worker_cpus(two_sockets,2).empty());require(select_worker_cpus(two_sockets,4).empty());
 require(select_worker_cpus({{20,0,0},{22,0,1}},2).empty());
 auto invalid=two_sockets;invalid.back().core=-1;require(select_worker_cpus(invalid,3).empty());
 auto reversed=two_sockets;std::reverse(reversed.begin(),reversed.end());require(select_worker_cpus(reversed,3)==std::vector<int>({20,22,47}));
#if defined(__linux__)
 cpu_set_t before,after;CPU_ZERO(&before);CPU_ZERO(&after);require(sched_getaffinity(0,sizeof(before),&before)==0);
 int selected=-1;for(int cpu=0;cpu<CPU_SETSIZE;++cpu)if(CPU_ISSET(cpu,&before)){selected=cpu;break;}require(selected>=0);
 bool pinned=false,exact=false,invalid_unchanged=false;
 std::thread worker([&]{pinned=place_worker(selected);cpu_set_t actual;CPU_ZERO(&actual);exact=sched_getaffinity(0,sizeof(actual),&actual)==0&&CPU_COUNT(&actual)==1&&CPU_ISSET(selected,&actual);invalid_unchanged=!place_worker(-1)&&!place_worker(CPU_SETSIZE);});worker.join();
 require(pinned&&exact&&invalid_unchanged);require(sched_getaffinity(0,sizeof(after),&after)==0&&CPU_EQUAL(&before,&after));
#endif
 std::cout<<"WORKER_PLACEMENT_PASS topology selection inherited-mask parent-unchanged fallback\n";
}
