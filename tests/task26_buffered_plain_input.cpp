#include "../src/merge_r8qkeep/buffered_plain_input.hpp"
#include <array>
#include <fstream>
#include <iostream>
#include <zlib.h>
using namespace fqdup::r8qkeep;
struct Batch {struct Row {std::size_t offset,header,sequence,quality;};std::string bytes;std::vector<Row> rows;};
static std::vector<std::array<std::string,3>> reference(const std::string& path){
    gzFile file=gzopen(path.c_str(),"rb");if(!file)throw std::runtime_error("reference open");
    const auto line=[&](std::string& out){out.clear();char buffer[8192];for(;;){if(!gzgets(file,buffer,sizeof(buffer)))return !out.empty();const auto n=std::strlen(buffer);if(n&&buffer[n-1]=='\n'){out.append(buffer,n-1);return true;}out.append(buffer,n);}};
    std::vector<std::array<std::string,3>> out;std::array<std::string,3> row;std::string plus;
    try{while(line(row[0])){if(!line(row[1])||!line(plus)||!line(row[2]))throw std::runtime_error("truncated reference");out.push_back(row);}}catch(...){gzclose(file);throw;}gzclose(file);return out;
}
int main(){
    char name[]="/tmp/task26-plain-reader-XXXXXX";const char* directory=::mkdtemp(name);if(!directory)throw std::runtime_error("fixture directory");
    unsigned serial=0;
    const auto check_one=[&](const std::string& text){const auto path=std::string(directory)+"/case-"+std::to_string(serial++);{std::ofstream out(path,std::ios::binary);out.write(text.data(),text.size());}
        bool bad=false;std::vector<std::array<std::string,3>> expected;try{expected=reference(path);}catch(const std::runtime_error&){bad=true;}
        for(std::size_t limit:{1,13,4096}){auto input=BufferedPlainInput::open(path);if(!input)throw std::runtime_error("plain input open");bool failed=false;std::vector<std::array<std::string,3>> actual;
            try{for(;;){Batch batch;input->read(batch,limit);if(batch.rows.empty())break;for(const auto& row:batch.rows){const auto p=row.offset;actual.push_back({batch.bytes.substr(p,row.header),batch.bytes.substr(p+row.header,row.sequence),batch.bytes.substr(p+row.header+row.sequence,row.quality)});}}}catch(const std::runtime_error&){failed=true;}
            if(failed!=bad||(!bad&&actual!=expected))throw std::runtime_error("plain reader parity");}
    };
    const std::string prefix="@before/1\nACGT\n+\nIIII\n@before/2\nTGCA\n+\nIIII\n";
    const auto check=[&](const std::string& text){check_one(text);check_one(prefix+text);};
    for(unsigned field=0;field<4;++field)for(unsigned pos=0;pos<65;++pos){
        std::array<std::string,4> lines{{"@"+std::string(64,'h'),std::string(65,'A'),"+"+std::string(64,'p'),std::string(65,'I')}};
        lines[field][pos]=0;std::string text;for(const auto& line:lines)text+=line+"\n";check(text);
    }
    check("");check("@x\nAC\n+\nII\n");check("@x\nAC\n+\nII");check("@x\n\n+\n\n");check("@x\r\nAC\r\n+\r\nII\r\n");
    for(const auto& text:{"@x","@x\n","@x\nAC","@x\nAC\n+","@x\nAC\n+\n"})check(text);
    for(unsigned length:{8190,8191,8192,16383,1048575,1048576,1048577}){std::string text="@x\n"+std::string(length,'A')+"\n+\n"+std::string(length,'I')+"\n";check(text);
        for(auto at:{1U,length/2,length}){auto binary=text;binary[at]=0;check(binary);}}
    std::string text;for(unsigned i=0;i<12000;++i)text+="@"+std::to_string(i)+" comment\n"+std::string(i%512,'A')+"\n+ comment\n"+std::string(i%512,'I')+"\n";check(text);
    const auto compressed=std::string(directory)+"/gzip";gzFile file=gzopen(compressed.c_str(),"wb");gzputs(file,"@x\nAC\n+\nII\n");gzclose(file);
    if(BufferedPlainInput::open(compressed)||BufferedPlainInput::open("-")||BufferedPlainInput::open(directory))throw std::runtime_error("plain fallback");
    std::cout<<"buffered plain input parity: boundaries, NUL, CR, EOF, empty, truncation and gzip fallback PASS\n";
}
