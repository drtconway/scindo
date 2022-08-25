#ifndef SCINDO_FACTOR_HPP
#define SCINDO_FACTOR_HPP

namespace scindo {

struct level : std::tuple<size_t, std::string_view> {
    using std::tuple<size_t, std::string_view>::std::tuple<size_t, std::string_view>;
};

struct factor {
  const std::vector<std::string> levels;
  std::unordered_map<std::string_view, size_t> level_index;

  factor(const std::vector<std::string> &p_levels) : levels(p_levels) {
    for (size_t i = 0; i < levels.size(); ++i) {
        level_index[std::string_view(levels[i])] = i;
    }
  }

  factor(const factor&) = delete;

  level operator[](const std::string& p_label) const {
    auto itr = level_index.find(p_label);
    if (itr == level_index.end()) {
        throw  std::runtime_error("unable to find label in factor.");
    }
    return level(itr->second, itr->first);
  }
}
} // namespace scindo
// namespace scindo

#endif // SCINDO_FACTOR_HPP