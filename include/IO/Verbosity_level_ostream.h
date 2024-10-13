#ifndef IO_VERBOSITY_LEVEL_OSTREAM_H
#define IO_VERBOSITY_LEVEL_OSTREAM_H

#include <iostream>

namespace CGAL::IO {
class Verbosity_level_ostream {
   public:
    Verbosity_level_ostream(std::ostream &os = std::clog, int output_level = 0, int level = 0)
        : _os(os), _level(level), _output_level(output_level) {}

    int level() const { return _level; }

    int output_level() const { return _output_level; }

    /**
     * @brief Change the verbosity level of current stream.
     *
     * @param level The verbosity level.
     * @return Verbosity_level_ostream&
     */
    Verbosity_level_ostream &level(int level) {
        _level = level;
        return *this;
    }

    /**
     * @brief Set the maximum output level of verbosity. Any output with a level higher than this will be suppressed.
     *
     * @param output_level The maximum output level of verbosity.
     * @return Verbosity_level_ostream&
     */
    Verbosity_level_ostream &output_level(int output_level) {
        _output_level = output_level;
        return *this;
    }

    Verbosity_level_ostream &operator<<(Verbosity_level_ostream &(*f)(Verbosity_level_ostream &)) { return f(*this); }

    template <typename T>
    Verbosity_level_ostream &operator<<(T &&t) {
        if (_level <= _output_level) {
            _os << std::forward<T>(t);
        }
        return *this;
    }

    Verbosity_level_ostream &operator<<(std::ostream &(*f)(std::ostream &)) {
        if (_level <= _output_level) {
            f(_os);
        }
        return *this;
    }

    Verbosity_level_ostream &operator<<(std::ios &(*f)(std::ios &)) {
        if (_level <= _output_level) {
            f(_os);
        }
        return *this;
    }

   private:
    std::ostream &_os;
    int _level, _output_level;
};

struct Set_level {
    int level;
    Set_level(int level) : level(level) {}
};

inline Set_level level(int level) {
    return Set_level(level);
}

inline Verbosity_level_ostream &operator<<(Verbosity_level_ostream &l, Set_level sl) {
    l.level(sl.level);
    return l;
}

}  // namespace CGAL::IO

#endif  // IO_VERBOSITY_LEVEL_OSTREAM_H
