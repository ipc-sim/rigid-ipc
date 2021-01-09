#pragma once

namespace ccd {

template <typename T> inline auto get_type_name()
{
    if (std::is_same<T, Interval>::value) {
        return "Interval";
    } else if (std::is_same<T, double>::value) {
        return "double";
    }
    return typeid(T).name();
}

} // namespace ccd
