#pragma once
#include <algorithm>
#include <cstdlib>
#include <limits>
#include <thread>
#include <fstream>
#include <map>
#include <string>
#include <utility>
#include <vector>
#if defined(__linux__)
#include <sched.h>
#endif
namespace fqdup::r8qkeep {
// Cache the inherited allocation before any worker narrows its affinity.
inline int allocation_cpu_count(){
    static const int count=[] {
        int available=static_cast<int>(std::max(1u,std::thread::hardware_concurrency()));
#if defined(__linux__)
        cpu_set_t mask;CPU_ZERO(&mask);
        if(sched_getaffinity(0,sizeof(mask),&mask)==0)available=std::max(1,CPU_COUNT(&mask));
#endif
        if(const char* value=std::getenv("SLURM_CPUS_PER_TASK")){
            char* end=nullptr;const long allocated=std::strtol(value,&end,10);
            if(end!=value&&*end=='\0'&&allocated>0&&allocated<=std::numeric_limits<int>::max())
                available=std::min(available,static_cast<int>(allocated));
        }
        return available;
    }();
    return count;
}
inline int bounded_worker_count(int requested){return std::max(1,std::min(requested,allocation_cpu_count()));}
struct WorkerCpu {int logical,package,core;};
// Only specialize a CPU set that supplies exactly one physical core per
// worker plus SMT siblings. Larger, smaller and incomplete sets stay unbound.
inline std::vector<int> select_worker_cpus(const std::vector<WorkerCpu>& allowed,int workers){
    if(workers<2||allowed.size()<=std::size_t(workers))return {};
    std::map<std::pair<int,int>,int> first;
    for(const auto cpu:allowed){
        if(cpu.logical<0||cpu.package<0||cpu.core<0)return {};
        const auto key=std::make_pair(cpu.package,cpu.core);
        auto inserted=first.emplace(key,cpu.logical);
        if(!inserted.second)inserted.first->second=std::min(inserted.first->second,cpu.logical);
    }
    if(first.size()!=std::size_t(workers))return {};
    std::vector<int> result;for(const auto& entry:first)result.push_back(entry.second);
    std::sort(result.begin(),result.end());return result;
}
inline std::vector<int> worker_cpu_plan(int workers){
    if(workers<2)return {};
#if defined(__linux__)
    cpu_set_t mask;CPU_ZERO(&mask);
    if(sched_getaffinity(0,sizeof(mask),&mask)!=0||CPU_COUNT(&mask)<=workers)return {};
    std::vector<WorkerCpu> allowed;
    for(int cpu=0;cpu<CPU_SETSIZE;++cpu)if(CPU_ISSET(cpu,&mask)){
        const auto base="/sys/devices/system/cpu/cpu"+std::to_string(cpu)+"/topology/";
        int package=-1,core=-1;
        std::ifstream package_file(base+"physical_package_id"),core_file(base+"core_id");
        if(!(package_file>>package)||!(core_file>>core))return {};
        allowed.push_back({cpu,package,core});
    }
    return select_worker_cpus(allowed,workers);
#else
    (void)workers;return {};
#endif
}
inline bool place_worker(int cpu){
#if defined(__linux__)
    if(cpu<0||cpu>=CPU_SETSIZE)return false;
    cpu_set_t inherited;CPU_ZERO(&inherited);
    if(sched_getaffinity(0,sizeof(inherited),&inherited)!=0||!CPU_ISSET(cpu,&inherited))return false;
    cpu_set_t target;CPU_ZERO(&target);CPU_SET(cpu,&target);
    return sched_setaffinity(0,sizeof(target),&target)==0;
#else
    (void)cpu;return false;
#endif
}
}
