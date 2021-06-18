#ifndef SCINDO_PROFILE_HPP
#define SCINDO_PROFILE_HPP

#include <chrono>
#include <memory>
#include <unordered_map>
#include <boost/log/trivial.hpp>

namespace scindo
{
    template <bool Enabled = true>
    struct profile
    {
        struct profile_record
        {
            const char* label;
            size_t num_calls;
            std::chrono::duration<double> time;

            profile_record(const char* p_label)
                : label(p_label), num_calls(0), time(0)
            {
            }
        };
        using profile_record_ptr = std::shared_ptr<profile_record>;

        profile_record_ptr rec;
        std::chrono::time_point<std::chrono::high_resolution_clock> start;

        explicit profile(const char* p_label)
        {
            if (Enabled)
            {
                rec = get_rec(p_label);
                start = std::chrono::high_resolution_clock::now();
            }
        }

        ~profile()
        {
            if (Enabled)
            {
                auto finish = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> delta = (finish - start);
                rec->num_calls += 1;
                rec->time += delta;
                if (as_you_go())
                {
                    BOOST_LOG_TRIVIAL(info)
                        << rec->label << '\t' << rec->num_calls << '\t' << delta.count() << '\t' << rec->time.count();
                }
            }
        }

        static profile_record_ptr get_rec(const char* p_label)
        {
            profile_record_ptr& ptr = records()[p_label];
            if (!ptr)
            {
                ptr = profile_record_ptr(new profile_record(p_label));
            }
            return ptr;
        }

        static std::unordered_map<const char*, profile_record_ptr>& records()
        {
            static std::unordered_map<const char*, profile_record_ptr> R;
            return R;
        }

        static bool& as_you_go()
        {
            static bool a = false;
            return a;
        }

        static void report()
        {
            for (auto itr = records().begin(); itr != records().end(); ++itr)
            {
                const auto& rec = itr->second;
                BOOST_LOG_TRIVIAL(info)
                    << rec->label << '\t' << rec->num_calls << '\t' << rec->time.count() << '\t' << (rec->time.count() / rec->num_calls);
            }
        }
    };
}
// namespace scindo

#endif // SCINDO_PROFILE_HPP
