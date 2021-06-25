#ifndef SCINDO_BAM_HPP
#define SCINDO_BAM_HPP

#include <stdexcept>

extern "C"
{
#include <htslib/hts.h>
#include <htslib/sam.h>
}

namespace scindo
{
    namespace bam
    {
        struct seq
        {
            size_t z;
            const uint8_t* ptr;
            size_t off;

            seq(size_t p_z, const uint8_t* p_ptr, size_t p_offset = 0)
                : z(p_z), ptr(p_ptr), off(p_offset)
            {
            }

            size_t size() const
            {
                return z;
            }

            char operator[](size_t p_idx) const
            {
                uint8_t x = bam_seqi(ptr, off + p_idx);
                return "*AC3G567T9abcdeN"[x];
            }

            seq range(size_t p_begin, size_t p_end) const
            {
                return seq(p_end - p_begin, ptr, off + p_begin);
            }

            std::string asString() const
            {
                std::string res;
                res.reserve(size());
                for (size_t i = 0; i < size(); ++i)
                {
                    res.push_back((*this)[i]);
                }
                return res;
            }
        };

        struct qual
        {
            size_t z;
            const char* ptr;

            qual(size_t p_z, const char* p_ptr)
                : z(p_z), ptr(p_ptr)
            {
            }

            size_t size() const
            {
                return z;
            }

            char operator[](size_t p_idx) const
            {
                return ptr[p_idx];
            }
        };

        struct cigar
        {
            size_t z;
            const uint32_t* ptr;

            cigar(size_t p_z, const uint32_t* p_ptr)
                : z(p_z), ptr(p_ptr)
            {
            }

            size_t size() const
            {
                return z;
            }

            char op(size_t p_idx) const
            {
                return bam_cigar_opchr(ptr[p_idx]);
            }

            size_t len(size_t p_idx) const
            {
                return bam_cigar_oplen(ptr[p_idx]);
            }
        };

        struct flag
        {
            static constexpr uint16_t paired = BAM_FPAIRED;
            static constexpr uint16_t proper_pair = BAM_FPROPER_PAIR;;
            static constexpr uint16_t unmapped = BAM_FUNMAP;
            static constexpr uint16_t mate_unmapped = BAM_FMUNMAP;
            static constexpr uint16_t reverse = BAM_FREVERSE;
            static constexpr uint16_t mate_reverse = BAM_FMREVERSE;
            static constexpr uint16_t read1 = BAM_FREAD1;
            static constexpr uint16_t read2 = BAM_FREAD2;
            static constexpr uint16_t secondary = BAM_FSECONDARY;
            static constexpr uint16_t qcfail = BAM_FQCFAIL;
            static constexpr uint16_t duplicate = BAM_FDUP;
            static constexpr uint16_t supplementary = BAM_FSUPPLEMENTARY;

            uint16_t value;

            flag(uint16_t p_value)
                : value(p_value)
            {
            }

            template<uint16_t Which>
            bool is() const
            {
                return value & Which;
            }
        };

        struct bam_file_reader
        {
            const std::string filename;

            bam_file_reader(const std::string& p_filename)
                : filename(p_filename), file(sam_open(p_filename.c_str(), "r")),
                  header(NULL), aln(NULL)
            {
                if (!file)
                {
                    auto msg = std::string("unable to open BAM file `") + filename + std::string("` for reading");
                    throw std::runtime_error(msg);
                }

                header = sam_hdr_read(file);
                if (!header)
                {
                    auto msg = std::string("unable read header for BAM file `") + filename + std::string("`");
                    throw std::runtime_error(msg);
                }

                aln = bam_init1();
                if (!aln)
                {
                    auto msg = std::string("unable create alignment for BAM file ") + filename + std::string("`");
                    throw std::runtime_error(msg);
                }
            }

            ~bam_file_reader()
            {
                if (file)
                {
                    sam_close(file);
                    file = NULL;
                }
                if (header)
                {
                    sam_hdr_destroy(header);
                    header = NULL;
                }
                if (aln)
                {
                    bam_destroy1(aln);
                }
            }

            bool next()
            {
                int res = bam_read1(file->fp.bgzf, aln);
                if (res < -1)
                {
                    auto msg = std::string("error reading BAM file `") + filename + std::string("`");
                    throw std::runtime_error(msg);
                }
                return (res >= 0);
            }

            uint16_t flag() const
            {
                return aln->core.flag;
            }

            const std::string& chrom()
            {
                if (!chroms.contains(aln->core.tid))
                {
                    chroms[aln->core.tid] = std::string(sam_hdr_tid2name(header, aln->core.tid));
                }
                return chroms.at(aln->core.tid);
            }

            int64_t pos0() const
            {
                return aln->core.pos;
            }

            int64_t pos1() const
            {
                return aln->core.pos + 1;
            }

            int mapq() const
            {
                return aln->core.qual;
            }

            scindo::bam::cigar cigar() const
            {
                const uint32_t* p = bam_get_cigar(aln);
                size_t l = aln->core.n_cigar;
                return scindo::bam::cigar(l, p);
            }

            scindo::bam::seq seq() const
            {
                const uint8_t* p = bam_get_seq(aln);
                size_t l = aln->core.l_qseq;
                return scindo::bam::seq(l, p);
            }

            scindo::bam::qual qual() const
            {
                const char* p = reinterpret_cast<const char*>(bam_get_qual(aln));
                size_t l = aln->core.l_qseq;
                return scindo::bam::qual(l, p);
            }

            samFile* file;
            sam_hdr_t* header;
            bam1_t* aln;
            std::unordered_map<int,std::string> chroms;
        };

    }
    // namespace bam
}
// namespace scindo

#endif // SCINDO_BAM_HPP
