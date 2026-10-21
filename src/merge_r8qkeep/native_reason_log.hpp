#pragma once
#ifdef FQDUP_NATIVE_REASONS
#include <mutex>
namespace fqdup::r8qkeep {
inline std::mutex& native_reason_mutex() {
    static std::mutex mutex;
    return mutex;
}
}
#endif
