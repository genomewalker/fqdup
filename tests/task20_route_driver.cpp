#include "../src/merge_r8qkeep/repair.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
std::string read(const char* path){std::ifstream f(path);if(!f)throw std::runtime_error("input open");return {std::istreambuf_iterator<char>(f),{}};}
int main(int argc,char** argv){try{if(argc!=5)throw std::runtime_error("usage: route_driver NOMINAL_PRIOR FULL_PRIOR NATIVE.tsv NEW_OUTPUT_DIR");if(!std::filesystem::create_directory(argv[4]))throw std::runtime_error("output directory exists");fqdup::r8qkeep::RepairPipeline pipeline(read(argv[1]),read(argv[2]));const auto result=pipeline.route(read(argv[3]));const std::string out=argv[4];std::ofstream(out+"/routed.tsv")<<result.routed_tsv;std::ofstream(out+"/route-audit.tsv")<<result.audit_tsv;return 0;}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
