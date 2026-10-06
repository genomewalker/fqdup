#pragma once
#include <string>
namespace fqdup::r8qkeep {
// ECMAScript Number::toString decimal layout over shortest round-trip digits.
// Gate G3 remains authoritative for full prior-text parity.
std::string js_number(double value);
}
