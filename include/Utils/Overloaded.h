#ifndef UTILS_OVERLOADED_H
#define UTILS_OVERLOADED_H

template <class... Ts>
struct overloaded : Ts... {
    using Ts::operator()...;
};

// explicit deduction guide (not needed as of C++20)
template <class... Ts>
overloaded(Ts...) -> overloaded<Ts...>;

#endif  // UTILS_OVERLOADED_H
