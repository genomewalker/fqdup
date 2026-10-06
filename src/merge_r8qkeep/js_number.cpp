#include "js_number.hpp"
#include <charconv>
#include <cmath>
#include <stdexcept>
#include <string>
namespace fqdup::r8qkeep {
std::string js_number(double value) {
    if (std::isnan(value)) return "NaN";
    if (std::isinf(value)) return std::signbit(value) ? "-Infinity" : "Infinity";
    if (value == 0) return "0";
    const bool negative = value < 0;
    char buffer[128];
    const auto result = std::to_chars(buffer, buffer + sizeof(buffer), std::abs(value), std::chars_format::general);
    if (result.ec != std::errc()) throw std::runtime_error("shortest double conversion failed");
    std::string text(buffer, result.ptr);
    const auto e = text.find_first_of("eE");
    const int exponent = e == std::string::npos ? 0 : std::stoi(text.substr(e + 1));
    const auto mantissa = text.substr(0, e);
    const auto point = mantissa.find('.');
    int n = static_cast<int>(point == std::string::npos ? mantissa.size() : point) + exponent;
    std::string digits;
    for (char c : mantissa) if (c != '.') digits += c;
    while (digits.size() > 1 && digits.front() == '0') { digits.erase(0, 1); --n; }
    while (digits.size() > 1 && digits.back() == '0') digits.pop_back();
    const int k = static_cast<int>(digits.size());
    std::string out;
    if (k <= n && n <= 21) out = digits + std::string(n - k, '0');
    else if (0 < n && n <= 21) out = digits.substr(0, n) + "." + digits.substr(n);
    else if (-6 < n && n <= 0) out = "0." + std::string(-n, '0') + digits;
    else {
        out = digits.substr(0, 1);
        if (k > 1) out += "." + digits.substr(1);
        out += "e";
        if (n - 1 >= 0) out += "+";
        out += std::to_string(n - 1);
    }
    return negative ? "-" + out : out;
}
}
