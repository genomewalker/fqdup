#pragma once
// Shortest round-trip formatting with ECMAScript decimal layout keeps prior files byte-reproducible.
#include <string>
namespace fqdup::r8qkeep {
// ECMAScript Number::toString decimal layout over shortest round-trip digits.
// Gate G3 remains authoritative for full prior-text parity.
std::string shortest_double(double value);
}
