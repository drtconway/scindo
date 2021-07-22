#ifndef SCINDO_CONCURRENT_DEQUE_HPP
#define SCINDO_CONCURRENT_DEQUE_HPP

#include <deque>
#include <mutex>
#include <condition_variable>

namespace scindo
{
    template <typename T>
    class concurrent_deque
    {
    public:
        concurrent_deque(size_t p_max_size = 0)
            : m_max_size(p_max_size), m_ended(false)
        {
        }

        size_t size() const
        {
            return m_deque.size();
        }

        void end()
        {
            std::unique_lock<std::mutex> lk(m_lock);
            m_ended = true;
            m_consumer_cond.notify_all();
        }

        concurrent_deque& push_front(const T& p_item)
        {
            std::unique_lock<std::mutex> lk(m_lock);

            if (m_ended)
            {
                throw std::runtime_error("cannot push on ended concurrent_deque");
            }

            while (m_max_size > 0 && m_deque.size() >= m_max_size)
            {
                m_producer_cond.wait(lk);
            }

            m_deque.push_back(p_item);

            m_consumer_cond.notify_one();

            return *this;
        }

        concurrent_deque& push_front(T&& p_item)
        {
            std::unique_lock<std::mutex> lk(m_lock);

            if (m_ended)
            {
                throw std::runtime_error("cannot push on ended concurrent_deque");
            }

            while (m_max_size > 0 && m_deque.size() >= m_max_size)
            {
                m_producer_cond.wait(lk);
            }

            m_deque.push_back(p_item);

            m_consumer_cond.notify_one();

            return *this;
        }

        concurrent_deque& push_back(const T& p_item)
        {
            std::unique_lock<std::mutex> lk(m_lock);

            if (m_ended)
            {
                throw std::runtime_error("cannot push on ended concurrent_deque");
            }

            while (m_deque.size() >= m_max_size)
            {
                m_producer_cond.wait(lk);
            }

            m_deque.push_back(p_item);

            m_consumer_cond.notify_one();

            return *this;
        }

        concurrent_deque& push_back(T&& p_item)
        {
            std::unique_lock<std::mutex> lk(m_lock);

            if (m_ended)
            {
                throw std::runtime_error("cannot push on ended concurrent_deque");
            }

            while (m_deque.size() >= m_max_size)
            {
                m_producer_cond.wait(lk);
            }

            m_deque.push_back(p_item);

            m_consumer_cond.notify_one();

            return *this;
        }

        bool pop_front(T& p_item)
        {
            std::unique_lock<std::mutex> lk(m_lock);

            while (!m_ended && m_deque.size() == 0)
            {
                m_consumer_cond.wait(lk);
            }

            if (m_ended && m_deque.size() == 0)
            {
                return false;
            }

            std::swap(p_item, m_deque.front());
            m_deque.pop_front();

            m_producer_cond.notify_one();
            return true;
        }

        bool pop_back(T& p_item)
        {
            std::unique_lock<std::mutex> lk(m_lock);

            while (!m_ended && m_deque.size() == 0)
            {
                m_consumer_cond.wait(lk);
            }

            if (m_ended && m_deque.size() == 0)
            {
                return false;
            }

            std::swap(p_item, m_deque.back());
            m_deque.pop_back();

            m_producer_cond.notify_one();
            return true;
        }

    private:
        const size_t m_max_size;
        bool m_ended;
        std::mutex m_lock;
        std::condition_variable m_producer_cond;
        std::condition_variable m_consumer_cond;
        std::deque<T> m_deque;
    };
}
// namespace scindo

#endif // SCINDO_CONCURRENT_DEQUE_HPP
