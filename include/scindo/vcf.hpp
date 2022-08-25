#ifndef SCINDO_VCF_HPP
#define SCINDO_VCF_HPP

#include <cmath>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

extern "C" {
#include <htslib/hts.h>
#include <htslib/vcf.h>
}

namespace scindo {
namespace vcf {
struct allele {
  const int32_t val;

  allele(int32_t p_val) : val(p_val) {}

  bool is_missing() const { return bcf_gt_is_missing(val); }

  bool is_phased() const { return bcf_gt_is_phased(val); }

  int32_t operator*() const { return bcf_gt_allele(val); }
};

struct vcf_file_reader {
  struct info_id {
    int32_t id;
  };

  struct format_id {
    int32_t id;
  };

  const std::string filename;

  vcf_file_reader(const std::string &p_filename)
      : filename(p_filename), file(vcf_open(p_filename.c_str(), "r")),
        header(NULL), var(NULL) {
    if (!file) {
      auto msg = std::string("unable to open VCF file `") + filename +
                 std::string("` for reading");
      throw std::runtime_error(msg);
    }

    header = bcf_hdr_read(file);
    if (!header) {
      auto msg = std::string("unable read header for VCF file `") + filename +
                 std::string("`");
      throw std::runtime_error(msg);
    }
  }

  ~vcf_file_reader() {
    if (file) {
      vcf_close(file);
      file = NULL;
    }
    if (header) {
      bcf_hdr_destroy(header);
      header = NULL;
    }
    if (var) {
      bcf_destroy(var);
      var = NULL;
    }
  }

  std::vector<std::string> samples() const {
    std::vector<std::string> res;
    const size_t n = bcf_hdr_nsamples(header);
    for (size_t i = 0; i < n; ++i) {
      res.push_back(std::string(header->samples[i]));
    }
    return res;
  }

  bool lookup(const char *p_str, info_id &p_id) {
    p_id.id = bcf_hdr_id2int(header, BCF_DT_ID, p_str);
    return p_id.id >= 0;
  }

  bool lookup(const char *p_str, format_id &p_id) {
    p_id.id = bcf_hdr_id2int(header, BCF_DT_ID, p_str);
    return p_id.id >= 0;
  }

  bool next() {
    if (!var) {
      var = bcf_init();
      if (!var) {
        auto msg = std::string("unable to allocate variant for VCF file `") +
                   filename + std::string("`");
        throw std::runtime_error(msg);
      }
    }
    int res = bcf_read(file, header, var);
    if (res < -1) {
      auto msg =
          std::string("error reading VCF file `") + filename + std::string("`");
      throw std::runtime_error(msg);
    }
    return (res == 0);
  }

  const int chromId() const { return var->rid; }

  const std::string &chrom() {
    if (!chroms.contains(var->rid)) {
      chroms[var->rid] = std::string(bcf_hdr_id2name(header, var->rid));
    }
    return chroms.at(var->rid);
  }

  int64_t pos0() const { return var->pos; }

  int64_t pos1() const { return var->pos + 1; }

  std::string id() {
    bcf_unpack(var, BCF_UN_STR);
    return std::string(var->d.id);
  }

  std::string ref() {
    bcf_unpack(var, BCF_UN_STR);
    return std::string(var->d.allele[0]);
  }

  size_t num_alts() {
    bcf_unpack(var, BCF_UN_STR);
    return var->n_allele - 1;
  }

  std::string alt(size_t p_n) {
    bcf_unpack(var, BCF_UN_STR);
    return std::string(var->d.allele[p_n]);
  }

  float qual() const { return var->qual; }

  size_t num_filters() {
    bcf_unpack(var, BCF_UN_FLT);
    return var->d.n_flt;
  }

  const std::string &filter(size_t p_n) {
    bcf_unpack(var, BCF_UN_FLT);
    int32_t ix = var->d.flt[p_n];
    if (!filters.contains(ix)) {
      filters[ix] = std::string(header->id[BCF_DT_ID][ix].key);
    }
    return filters.at(ix);
  }

  size_t num_infos() {
    bcf_unpack(var, BCF_UN_INFO);
    return var->n_info;
  }

  const std::string &info_key(size_t p_n) {
    bcf_unpack(var, BCF_UN_INFO);
    const bcf_info_t &ifo = var->d.info[p_n];
    int32_t ix = ifo.key;
    if (!info_names.contains(ix)) {
      info_names[ix] = std::string(header->id[BCF_DT_ID][ix].key);
    }
    return info_names.at(ix);
  }

  template <typename T>
  std::enable_if<std::is_integral<T>::value, bool>::type info_value(size_t p_n,
                                                                    T &p_res) {
    bcf_unpack(var, BCF_UN_INFO);
    const bcf_info_t &ifo = var->d.info[p_n];
    p_res = ifo.v1.i;
    switch (ifo.type) {
    case BCF_BT_INT8:
    case BCF_BT_INT16:
    case BCF_BT_INT32:
    case BCF_BT_INT64:
    case BCF_BT_CHAR: {
      return true;
    }
    default: {
      return false;
    }
    }
  }

  template <typename T>
  std::enable_if<std::is_floating_point<T>::value, bool>::type
  info_value(size_t p_n, T &p_res) {
    bcf_unpack(var, BCF_UN_INFO);
    const bcf_info_t &ifo = var->d.info[p_n];
    p_res = ifo.v1.f;
    return ifo.type == BCF_BT_FLOAT;
  }

  bool info_value(size_t p_n, std::string &p_res) {
    bcf_unpack(var, BCF_UN_INFO);
    const bcf_info_t &ifo = var->d.info[p_n];
    const char *s = reinterpret_cast<const char *>(ifo.vptr);
    p_res.clear();
    p_res.insert(p_res.end(), s, s + ifo.len);
    return ifo.type == BCF_BT_NULL;
  }

  template <typename T> T info_value(size_t p_n) {
    T x;
    info_value(p_n, x);
    return x;
  }

  template <typename T>
  std::enable_if<std::is_integral<T>::value, bool>::type
  info_value(const info_id &p_id, T &p_res) {
    bcf_unpack(var, BCF_UN_INFO);
    const bcf_info_t *ptr = bcf_get_info_id(var, p_id.id);
    if (!ptr) {
      return false;
    }
    const bcf_info_t &ifo = *ptr;
    p_res = ifo.v1.i;
    switch (ifo.type) {
    case BCF_BT_INT8:
    case BCF_BT_INT16:
    case BCF_BT_INT32:
    case BCF_BT_INT64:
    case BCF_BT_CHAR: {
      return true;
    }
    default: {
      return false;
    }
    }
  }

  template <typename T>
  std::enable_if<std::is_floating_point<T>::value, bool>::type
  info_value(const info_id &p_id, T &p_res) {
    bcf_unpack(var, BCF_UN_INFO);
    const bcf_info_t *ptr = bcf_get_info_id(var, p_id.id);
    if (!ptr) {
      return false;
    }
    const bcf_info_t &ifo = *ptr;
    p_res = ifo.v1.f;
    return ifo.type == BCF_BT_FLOAT;
  }

  bool info_value(const info_id &p_id, std::string &p_res) {
    bcf_unpack(var, BCF_UN_INFO);
    const bcf_info_t *ptr = bcf_get_info_id(var, p_id.id);
    if (!ptr) {
      return false;
    }
    const bcf_info_t &ifo = *ptr;
    const char *s = reinterpret_cast<const char *>(ifo.vptr);
    p_res.clear();
    p_res.insert(p_res.end(), s, s + ifo.len);
    return ifo.type == BCF_BT_NULL;
  }

  bool info_value(const info_id &p_id, bool &p_res) {
    bcf_unpack(var, BCF_UN_INFO);
    const bcf_info_t *ptr = bcf_get_info_id(var, p_id.id);
    p_res = (ptr != NULL);
    return p_res;
  }

  template <typename T> std::optional<T> info_value(const info_id &p_id) {
    T x;
    if (info_value(p_id, x)) {
      return x;
    }
    return {};
  }

  size_t num_formats() {
    bcf_unpack(var, BCF_UN_FMT);
    return var->n_fmt;
  }

  const std::string &format_key(size_t p_n) {
    bcf_unpack(var, BCF_UN_FMT);
    const bcf_fmt_t &fmt = var->d.fmt[p_n];
    int32_t ix = fmt.id;
    if (!format_names.contains(ix)) {
      format_names[ix] = std::string(header->id[BCF_DT_ID][ix].key);
    }
    return format_names.at(ix);
  }

  template <typename T, typename E = void> struct format_type_tag {};

  template <typename T>
  struct format_type_tag<
      T, typename std::enable_if<std::is_same<T, int32_t>::value, void>::type> {
    static constexpr int type_tag = BCF_HT_INT;
    static constexpr int32_t missing_value = bcf_int32_missing;
  };

  template <typename T>
  struct format_type_tag<
      T, typename std::enable_if<std::is_same<T, int64_t>::value, void>::type> {
    static constexpr int type_tag = BCF_HT_LONG;
    static constexpr int64_t missing_value = bcf_int64_missing;
  };

  template <typename T>
  struct format_type_tag<
      T, typename std::enable_if<std::is_same<T, float>::value, void>::type> {
    static constexpr int type_tag = BCF_HT_REAL;
    static constexpr float missing_value = NAN;
  };

  template <typename T>
  struct format_type_tag<
      T, typename std::enable_if<std::is_same<T, std::string>::value,
                                 void>::type> {
    static constexpr int type_tag = BCF_HT_STR;
    static constexpr const char *const missing_value = "";
  };

  template <typename T> struct format_data {
    static constexpr T missing_value = format_type_tag<T>::missing_value;

    int32_t arr_size;
    void *arr_ptr;

    format_data() : arr_size(0), arr_ptr(NULL) {}

    ~format_data() {
      if (arr_ptr) {
        free(arr_ptr);
      }
    }

    size_t size() const { return arr_size; }

    const T *data() const { return reinterpret_cast<const T *>(arr_ptr); }

    const T &operator[](size_t p_idx) const { return data()[p_idx]; }
  };

  template <typename T>
  bool format_values(const char *p_tag, format_data<T> &p_res) {
    int r0 =
        bcf_get_format_values(header, var, p_tag, &(p_res.arr_ptr),
                              &(p_res.arr_size), format_type_tag<T>::type_tag);
    return r0 >= 0;
  }

  std::string to_string() {
    kstring_t s;
    ks_initialize(&s);
    vcf_format(header, var, &s);
    std::string res(ks_c_str(&s));
    ks_free(&s);
    return res;
  }

private:
  htsFile *file;
  bcf_hdr_t *header;
  bcf1_t *var;
  std::unordered_map<int, std::string> chroms;
  std::unordered_map<int32_t, std::string> filters;
  std::unordered_map<int32_t, std::string> info_names;
  std::unordered_map<int32_t, std::string> format_names;
};
} // namespace vcf
  // namespace vcf
} // namespace scindo
// namespace scindo

#endif // SCINDO_VCF_HPP
