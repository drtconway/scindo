#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "scindo/vcf.hpp"

#include <filesystem>

namespace // anonymous
{

}
// namespace anonymous

TEST_CASE("test template instantiation", "[vcf]") {
    REQUIRE( scindo::vcf::vcf_file_reader::format_type_tag<int32_t>::type_tag == BCF_HT_INT );
    REQUIRE( scindo::vcf::vcf_file_reader::format_type_tag<int64_t>::type_tag == BCF_HT_LONG );
    REQUIRE( scindo::vcf::vcf_file_reader::format_type_tag<float>::type_tag == BCF_HT_REAL );
    REQUIRE( scindo::vcf::vcf_file_reader::format_type_tag<std::string>::type_tag == BCF_HT_STR );
}

TEST_CASE("parse example 1 from the VCF spec", "[vcf]") {
    scindo::vcf::vcf_file_reader V("data/tiny-vcf.vcf");
    SECTION("samples") {
        std::vector<std::string> s = V.samples();
        REQUIRE( s.size() == 3 );
    }
    scindo::vcf::vcf_file_reader::info_id db;
    scindo::vcf::vcf_file_reader::info_id h2;
    SECTION("info id") {
        REQUIRE( V.lookup("DB", db) );
        REQUIRE( V.lookup("H2", h2) );
    }
    SECTION("variants") {
        size_t n = 0;
        while (V.next())
        {
            switch (n)
            {
                case 0:
                {
                    REQUIRE( V.chrom() == "20" );
                    REQUIRE( V.pos1() == 14370 );
                    REQUIRE( V.id() == "rs6054257" );
                    REQUIRE( V.ref() == "G" );
                    REQUIRE( V.num_alts() == 1 );
                    REQUIRE( V.alt(1) == "A" );
                    REQUIRE( V.qual() == 29 );
                    REQUIRE( V.num_filters() == 1 );
                    REQUIRE( V.filter(0) == "PASS" );
                    REQUIRE( V.num_infos() == 5 );
                    REQUIRE( V.info_key(0) == "NS" );
                    REQUIRE( V.info_value<int>(0) == 3 );
                    REQUIRE( V.info_key(1) == "DP" );
                    REQUIRE( V.info_value<int>(1) == 14 );
                    REQUIRE( V.info_key(2) == "AF" );
                    REQUIRE( V.info_value<float>(2) == 0.5 );
                    REQUIRE( V.info_key(3) == "DB" );
                    REQUIRE( V.info_key(4) == "H2" );
                    break;
                }
                case 1:
                {
                    REQUIRE( V.chrom() == "20" );
                    REQUIRE( V.pos1() == 17330 );
                    REQUIRE( V.id() == "." );
                    REQUIRE( V.ref() == "T" );
                    REQUIRE( V.num_alts() == 1 );
                    REQUIRE( V.alt(1) == "A" );
                    REQUIRE( V.qual() == 3 );
                    REQUIRE( V.num_filters() == 1 );
                    REQUIRE( V.filter(0) == "q10" );
                    break;
                }
                case 2:
                {
                    REQUIRE( V.chrom() == "20" );
                    REQUIRE( V.pos1() == 1110696 );
                    REQUIRE( V.id() == "rs6040355" );
                    REQUIRE( V.ref() == "A" );
                    REQUIRE( V.num_alts() == 2 );
                    REQUIRE( V.alt(1) == "G" );
                    REQUIRE( V.alt(2) == "T" );
                    REQUIRE( V.qual() == 67 );
                    REQUIRE( V.num_filters() == 1 );
                    REQUIRE( V.filter(0) == "PASS" );
                    REQUIRE( (bool)V.info_value<bool>(db) );
                    REQUIRE( !(bool)V.info_value<bool>(h2) );

                    scindo::vcf::vcf_file_reader::format_data<int32_t> gt;
                    REQUIRE( V.format_values("GT", gt) );
                    REQUIRE( gt.size() == 6 );
                    std::vector<scindo::vcf::allele> A(gt.data(), gt.data() + gt.size());
                    REQUIRE( !A[0].is_missing() );
                    REQUIRE( !A[0].is_phased() );
                    REQUIRE( *A[0] == 1 );
                    REQUIRE( !A[1].is_missing() );
                    REQUIRE( A[1].is_phased() );
                    REQUIRE( *A[1] == 2 );
                    REQUIRE( !A[2].is_missing() );
                    REQUIRE( !A[2].is_phased() );
                    REQUIRE( *A[2] == 2 );
                    REQUIRE( !A[3].is_missing() );
                    REQUIRE( A[3].is_phased() );
                    REQUIRE( *A[3] == 1 );
                    REQUIRE( !A[4].is_missing() );
                    REQUIRE( !A[4].is_phased() );
                    REQUIRE( *A[4] == 2 );
                    REQUIRE( !A[5].is_missing() );
                    REQUIRE( !A[5].is_phased() );
                    REQUIRE( *A[5] == 2 );
                    break;
                }
                case 3:
                {
                    REQUIRE( V.chrom() == "20" );
                    REQUIRE( V.pos1() == 1230237 );
                    scindo::vcf::vcf_file_reader::format_data<int32_t> gt;
                    REQUIRE( V.format_values("GT", gt) );
                    REQUIRE( gt.size() == 6 );
                    std::vector<scindo::vcf::allele> A(gt.data(), gt.data() + gt.size());
                    scindo::vcf::vcf_file_reader::format_data<int32_t> hq;
                    REQUIRE( V.format_values("HQ", hq) );
                    REQUIRE( hq.size() == 6 );
                    REQUIRE( hq[0] == 56 );
                    REQUIRE( hq[1] == 60 );
                    REQUIRE( hq[2] == 51 );
                    REQUIRE( hq[3] == 51 );
                    REQUIRE( hq[4] == scindo::vcf::vcf_file_reader::format_data<int32_t>::missing_value );
                    break;
                }
                case 4:
                {
                    REQUIRE( V.chrom() == "20" );
                    REQUIRE( V.pos1() == 1234567 );
                    break;
                }
                default:
                {
                    REQUIRE( false );
                }
            }
            n += 1;
        }
        REQUIRE( n == 5 );
    }
}
