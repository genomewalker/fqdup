#pragma once
#include <array>
#include <sstream>
#include <string>
#include <vector>
namespace fqdup::r8qkeep {
inline constexpr std::array<double,34> cutoff_grid{{.0001,.000133,.000178,.000237,.000316,.000422,.000562,.00075,.001,.00133,.00178,.00237,.00316,.00422,.00562,.0075,.01,.0133,.0178,.0237,.0316,.0422,.05,.0562,.075,.1,.133,.178,.237,.316,.422,.562,.75,1}};
inline std::vector<std::string> cutoff_grid_labels(){std::vector<std::string> labels;for(double c:cutoff_grid){std::ostringstream text;text<<c;labels.push_back(text.str());}return labels;}
}
