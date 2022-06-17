#include "common.hpp"

ScopedTimer::ScopedTimer(const char *func, bool in_human_readable)
    : human_readable(in_human_readable), function_name_{ func }, start_{
          ClockType::now()
      }
{
}
ScopedTimer::~ScopedTimer()
{
    using namespace std::chrono;
    using days = duration<int, std::ratio<86400> >;

    auto stop = ClockType::now();
    auto duration = (stop - start_);
    auto ms = duration_cast<milliseconds>(duration);
    if (human_readable) {
        auto d = duration_cast<days>(duration);
        duration -= d;
        auto h = duration_cast<hours>(duration);
        duration -= h;
        auto m = duration_cast<minutes>(duration);
        duration -= m;
        auto s = duration_cast<seconds>(duration);
        duration -= s;
        auto ms = duration_cast<milliseconds>(duration);
        if (d.count() > 0)
            std::cerr << d.count() << "days ";
        if (d.count() || h.count())
            std::cerr << h.count() << "hrs ";
        if (d.count() || h.count() || m.count())
            std::cerr << m.count() << "min ";
        if (d.count() || h.count() || m.count() || s.count())
            std::cerr << s.count() << "sec ";
        std::cerr << ms.count() << "ms\t" << function_name_ << '\n';
    } else {
        std::cerr << ms.count() << "ms\t" << function_name_ << '\n';
    }
}
