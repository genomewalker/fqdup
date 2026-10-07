#pragma once
#include <streambuf>
#include <string>
namespace fqdup::r8qkeep {
// istream get area over an immutable, caller-owned string; no whole-text copy.
class StringInputBuffer final:public std::streambuf {
public:
    explicit StringInputBuffer(const std::string& text){auto* begin=const_cast<char*>(text.data());setg(begin,begin,begin+text.size());}
};
}
