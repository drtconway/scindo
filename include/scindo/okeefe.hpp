#ifndef SCINDO_OKEEFE_HPP
#define SCINDO_OKEEFE_HPP

#include <vector>

namespace scindo {

template <typename T> struct okeefe {
  std::vector<T> queue_front;
  std::vector<T> queue_back;
  size_t max_size;

  okeefe() : max_size(0) {}

  size_t size() const { return queue_front.size() + queue_back.size(); }

  void reserve(size_t p_n) { queue_back.reserve(p_n); }

  const T &front() const {
    if (queue_front.size() > 0) {
      return queue_front.back();
    } else {
      return queue_back.front();
    }
  }
  T &front() {
    if (queue_front.size() > 0) {
      return queue_front.back();
    } else {
      return queue_back.front();
    }
  }

  const T &back() const {
    if (queue_back.size() > 0) {
      return queue_back.back();
    } else {
      return queue_front.front();
    }
  }
  T &back() {
    if (queue_back.size() > 0) {
      return queue_back.back();
    } else {
      return queue_front.front();
    }
  }

  const T &operator[](size_t p_idx) const {
    if (p_idx < queue_front.size()) {
      return queue_front[queue_front.size() - 1 - p_idx];
    } else {
      p_idx -= queue_front.size();
      return queue_back[p_idx];
    }
  }
  T &operator[](size_t p_idx) {
    if (p_idx < queue_front.size()) {
      return queue_front[queue_front.size() - 1 - p_idx];
    } else {
      p_idx -= queue_front.size();
      return queue_back[p_idx];
    }
  }

  okeefe &push_front(const T &p_item) {
    queue_front.push_back(p_item);
    max_size = std::max(max_size, size());
    return *this;
  }

  okeefe &push_back(const T &p_item) {
    queue_back.push_back(p_item);
    max_size = std::max(max_size, size());
    return *this;
  }

  void pop_front() {
    ensure_front();
    queue_front.pop_back();
  }

  void pop_back() {
    ensure_back();
    queue_back.pop_back();
  }

  void ensure_front() {
    if (queue_front.size() == 0) {
      std::reverse(queue_back.begin(), queue_back.end());
      queue_back.swap(queue_front);
    }
  }

  void ensure_back() {
    if (queue_back.size() == 0) {
      std::reverse(queue_front.begin(), queue_front.end());
      queue_front.swap(queue_back);
    }
  }
};

} // namespace scindo

#endif // SCINDO_OKEEFE_HPP