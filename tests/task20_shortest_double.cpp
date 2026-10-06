#include "../src/merge_r8qkeep/shortest_double.hpp"
#include <cmath>
#include <iostream>
#include <limits>
#include <utility>
#include <vector>
int main() {
    using fqdup::r8qkeep::shortest_double;
    const std::vector<std::pair<double, const char*>> examples{
        {0.0,"0"},{-0.0,"0"},{1.0,"1"},{-1.25,"-1.25"},
        {0.316,"0.316"},{0.75,"0.75"},{1e-6,"0.000001"},{1e-7,"1e-7"},
        {1e20,"100000000000000000000"},{1e21,"1e+21"},{1.25e21,"1.25e+21"},
        {1e-15,"1e-15"},{1e-300,"1e-300"},
        {std::numeric_limits<double>::denorm_min(),"5e-324"},
        {std::numeric_limits<double>::max(),"1.7976931348623157e+308"},
        {std::nextafter(1.0,2.0),"1.0000000000000002"},
        {1000000000000000128.0,"1000000000000000100"}
    };
    for(const auto& [value,expected]:examples) {
        const auto actual=shortest_double(value);
        if(actual!=expected){std::cerr<<expected<<" != "<<actual<<'\n';return 1;}
        if(value!=std::numeric_limits<double>::denorm_min() && std::stod(actual)!=value)return 2;
    }
    if(shortest_double(std::numeric_limits<double>::quiet_NaN())!="NaN")return 3;
    if(shortest_double(std::numeric_limits<double>::infinity())!="Infinity")return 4;
    std::cout<<"JS_NUMBER_LAYOUT_PASS "<<examples.size()+2<<" fixtures\n";
}
