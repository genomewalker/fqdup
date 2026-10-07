#pragma once
#include <fstream>
#include <sstream>
#include <string>
namespace fqdup::r8qkeep {
inline std::string memory_status(){
    std::ifstream input("/proc/self/status");std::string line;long rss=-1,hwm=-1;
    while(std::getline(input,line)){std::istringstream fields(line);std::string key;long value;
        if(fields>>key>>value){if(key=="VmRSS:")rss=value;else if(key=="VmHWM:")hwm=value;}}
    return " rss_kib="+std::to_string(rss)+" hwm_kib="+std::to_string(hwm);
}
}
