#include "../src/merge_r8qkeep/streaming_output.hpp"
#include <fstream>
#include <iostream>
using namespace fqdup::r8qkeep;
static unsigned calls=0;static std::uint64_t scheduled=0;
static int record(int,std::uint64_t offset,std::uint64_t size){if(offset!=scheduled||!size)throw std::runtime_error("writeback range");scheduled+=size;++calls;return 0;}
static int unsupported(int,std::uint64_t,std::uint64_t){++calls;errno=EOPNOTSUPP;return -1;}
static int io_error(int,std::uint64_t,std::uint64_t){errno=EIO;return -1;}
int main(){char pattern[]="/tmp/fqdup-writeback-XXXXXX";const auto raw=::mkdtemp(pattern);if(!raw)throw std::runtime_error("temporary output");const std::filesystem::path temp(raw);
 struct Cleanup{std::filesystem::path path;~Cleanup(){std::error_code error;std::filesystem::remove_all(path,error);}}cleanup{temp};
 std::string expected(9*(1U<<20)+37,'\0');for(std::size_t i=0;i<expected.size();++i)expected[i]=char((i*31+i/251)%256);
 for(int mode=0;mode<3;++mode){calls=0;scheduled=0;StreamingOutput output=mode==0?StreamingOutput():StreamingOutput(65536,mode==1?record:unsupported);const auto path=temp/std::to_string(mode);output.open(path);if(!output)throw std::runtime_error("open");
  for(std::size_t at=0;at<expected.size();at+=40001)output.append(std::string_view(expected).substr(at,std::min<std::size_t>(40001,expected.size()-at)));
  StreamingOutput moved(std::move(output));moved.close();if(!moved)throw std::runtime_error("close");
  std::ifstream input(path,std::ios::binary);const std::string observed((std::istreambuf_iterator<char>(input)),{});if(expected!=observed)throw std::runtime_error("binary mismatch");if(mode==1&&calls<2)throw std::runtime_error("missing hints");if(mode==2&&calls!=1)throw std::runtime_error("unsupported fallback");
 }
 StreamingOutput failed(1,io_error);failed.open(temp/"error");bool caught=false;try{failed.append("x");}catch(const std::runtime_error&){caught=true;}failed.close();if(!caught||failed)throw std::runtime_error("writeback error not propagated");
 StreamingOutput missing;missing.open(temp/"missing"/"output");if(missing)throw std::runtime_error("open error ignored");
 std::cout<<"STREAMING_OUTPUT_PASS binary bytes ordered ranges unsupported fallback move ownership I/O errors\n";
}
