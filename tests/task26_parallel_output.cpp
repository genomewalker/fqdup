#include "../src/merge_r8qkeep/parallel_output.hpp"
#include <fstream>
#include <iostream>
#include <atomic>
using namespace fqdup::r8qkeep;
static std::promise<void> failed_write,allow_second;static auto second_ready=allow_second.get_future().share();static std::atomic<bool> invalid_fd{false};
static int fail(int,std::uint64_t,std::uint64_t){failed_write.set_value();errno=EIO;return -1;}
static int delayed(int fd,std::uint64_t,std::uint64_t){second_ready.wait();if(::fcntl(fd,F_GETFD)<0)invalid_fd=true;return 0;}
int main(){char pattern[]="/tmp/fqdup-parallel-output-XXXXXX";const auto raw=::mkdtemp(pattern);if(!raw)throw std::runtime_error("temporary output");const std::filesystem::path temp(raw);
 struct Cleanup{std::filesystem::path path;~Cleanup(){std::error_code error;std::filesystem::remove_all(path,error);}}cleanup{temp};
 std::array<StreamingOutput,5> files;std::array<std::string,5> expected;
 for(unsigned f=0;f<5;++f)files[f].open(temp/std::to_string(f));
 {
  ParallelOutputFiles writers(files,true,3);
  for(unsigned chunk=0;chunk<300;++chunk){std::array<std::string,5> data;std::array<std::string*,5> view;
   for(unsigned f=0;f<5;++f){if((chunk+f)%7)for(unsigned i=0,n=(chunk*913+f*37)%25001;i<n;++i)data[f]+=char((i+chunk*19+f*71)%256);expected[f]+=data[f];view[f]=&data[f];}
   writers.append(view);
  }
  writers.finish();if(writers.seconds()<0||writers.text_capacity_bound()==0)throw std::runtime_error("writer measurements");
 }
 for(unsigned f=0;f<5;++f){files[f].close();if(!files[f])throw std::runtime_error("close");std::ifstream input(temp/std::to_string(f),std::ios::binary);const std::string observed((std::istreambuf_iterator<char>(input)),{});if(observed!=expected[f])throw std::runtime_error("ordered bytes");}
 {
  std::array<StreamingOutput,5> outputs{StreamingOutput(1,fail),StreamingOutput(1,delayed),StreamingOutput(),StreamingOutput(),StreamingOutput()};
  for(unsigned f=0;f<5;++f)outputs[f].open(temp/("error-"+std::to_string(f)));
  auto release=std::thread([]{failed_write.get_future().wait();allow_second.set_value();});bool caught=false;
  {
   ParallelOutputFiles writers(outputs,false,2);std::array<std::string,5> text{"a","b","c","d","e"};writers.append({&text[0],&text[1],&text[2],&text[3],&text[4]});
   try{writers.finish();}catch(const std::runtime_error&){caught=true;}
  }
  release.join();for(auto& file:outputs)file.close();if(!caught||outputs[0]||invalid_fd)throw std::runtime_error("error drain ownership");
 }
 std::cout<<"PARALLEL_OUTPUT_PASS per-file binary order bounded queues empty chunks timing error drain descriptor lifetime\n";
}
