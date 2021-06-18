#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "scindo/gtf.hpp"

#include <sstream>

namespace // anonymous
{
}
// namespace anonymous

TEST_CASE("parse bits and pieces", "[gtf]") {
    SECTION("start") {
        const std::string l("2590639");
        tao::pegtl::memory_input in(l, "string");
        scindo::gtf::detail::state s;
        REQUIRE( tao::pegtl::parse<scindo::gtf::detail::integer, scindo::gtf::detail::action>(in, s) );
        REQUIRE( s.value.index() == 1 );
        REQUIRE( std::get<int64_t>(s.value) == 2590639 );
    }
    SECTION("attrs") {
        const std::string l("gene_id \"ENSG00000142606.16\"; gene_type \"protein_coding\"; gene_name \"MMEL1\"; level 2; hgnc_id \"HGNC:14668\"; havana_gene \"OTTHUMG00000000846.5\";");
        tao::pegtl::memory_input in(l, "string");
        scindo::gtf::detail::state s;
        REQUIRE( tao::pegtl::parse<scindo::gtf::detail::attrs, scindo::gtf::detail::action>(in, s) );
    }
    SECTION("seqname") {
        const std::string l("chr1");
        tao::pegtl::memory_input in(l, "string");
        scindo::gtf::detail::state s;
        REQUIRE( tao::pegtl::parse<scindo::gtf::detail::seqname, scindo::gtf::detail::action>(in, s) );
        REQUIRE( s.fields.contains(scindo::gtf::detail::state::seqname) );
        REQUIRE( s.fields.at(scindo::gtf::detail::state::seqname).index() == 0 );
        REQUIRE( std::get<std::string>(s.fields.at(scindo::gtf::detail::state::seqname)) == "chr1" );
    }
    SECTION("source") {
        const std::string l("HAVANA");
        tao::pegtl::memory_input in(l, "string");
        scindo::gtf::detail::state s;
        REQUIRE( tao::pegtl::parse<scindo::gtf::detail::source, scindo::gtf::detail::action>(in, s) );
        REQUIRE( s.fields.contains(scindo::gtf::detail::state::source) );
        REQUIRE( s.fields.at(scindo::gtf::detail::state::source).index() == 0 );
        REQUIRE( std::get<std::string>(s.fields.at(scindo::gtf::detail::state::source)) == "HAVANA" );
    }
    SECTION("feature") {
        const std::string l("gene");
        tao::pegtl::memory_input in(l, "string");
        scindo::gtf::detail::state s;
        REQUIRE( tao::pegtl::parse<scindo::gtf::detail::feature, scindo::gtf::detail::action>(in, s) );
        REQUIRE( s.fields.contains(scindo::gtf::detail::state::feature) );
        REQUIRE( s.fields.at(scindo::gtf::detail::state::feature).index() == 0 );
        REQUIRE( std::get<std::string>(s.fields.at(scindo::gtf::detail::state::feature)) == "gene" );
    }
    SECTION("start") {
        const std::string l("2590639");
        tao::pegtl::memory_input in(l, "string");
        scindo::gtf::detail::state s;
        REQUIRE( tao::pegtl::parse<scindo::gtf::detail::start, scindo::gtf::detail::action>(in, s) );
        REQUIRE( s.fields.contains(scindo::gtf::detail::state::start) );
        REQUIRE( s.fields.at(scindo::gtf::detail::state::start).index() == 1 );
        REQUIRE( std::get<int64_t>(s.fields.at(scindo::gtf::detail::state::start)) == 2590639 );
    }
    SECTION("end") {
        const std::string l("2633016");
        tao::pegtl::memory_input in(l, "string");
        scindo::gtf::detail::state s;
        REQUIRE( tao::pegtl::parse<scindo::gtf::detail::end, scindo::gtf::detail::action>(in, s) );
        REQUIRE( s.fields.contains(scindo::gtf::detail::state::end) );
        REQUIRE( s.fields.at(scindo::gtf::detail::state::end).index() == 1 );
        REQUIRE( std::get<int64_t>(s.fields.at(scindo::gtf::detail::state::end)) == 2633016 );
    }
    SECTION("score") {
        const std::string l("23.45");
        tao::pegtl::memory_input in(l, "string");
        scindo::gtf::detail::state s;
        REQUIRE( tao::pegtl::parse<scindo::gtf::detail::score, scindo::gtf::detail::action>(in, s) );
        REQUIRE( s.fields.contains(scindo::gtf::detail::state::score) );
        REQUIRE( s.fields.at(scindo::gtf::detail::state::score).index() == 2 );
        REQUIRE( std::get<double>(s.fields.at(scindo::gtf::detail::state::score)) == 23.45 );
    }
    SECTION("!score") {
        const std::string l(".");
        tao::pegtl::memory_input in(l, "string");
        scindo::gtf::detail::state s;
        REQUIRE( tao::pegtl::parse<tao::pegtl::sor<scindo::gtf::detail::score,scindo::gtf::detail::missing>, scindo::gtf::detail::action>(in, s) );
        REQUIRE( !s.fields.contains(scindo::gtf::detail::state::score) );
    }
    SECTION("strand") {
        const std::string l("-");
        tao::pegtl::memory_input in(l, "string");
        scindo::gtf::detail::state s;
        REQUIRE( tao::pegtl::parse<scindo::gtf::detail::strand, scindo::gtf::detail::action>(in, s) );
        REQUIRE( s.fields.contains(scindo::gtf::detail::state::strand) );
        REQUIRE( s.fields.at(scindo::gtf::detail::state::strand).index() == 0 );
        REQUIRE( std::get<std::string>(s.fields.at(scindo::gtf::detail::state::strand)) == "-" );
    }
    SECTION("frame") {
        const std::string l("0");
        tao::pegtl::memory_input in(l, "string");
        scindo::gtf::detail::state s;
        REQUIRE( tao::pegtl::parse<scindo::gtf::detail::frame, scindo::gtf::detail::action>(in, s) );
        REQUIRE( s.fields.contains(scindo::gtf::detail::state::frame) );
        REQUIRE( s.fields.at(scindo::gtf::detail::state::frame).index() == 0 );
        REQUIRE( std::get<std::string>(s.fields.at(scindo::gtf::detail::state::frame)) == "0" );
    }
    SECTION("full row") {
        const std::string l("chr1	HAVANA	gene	2590639	2633016	.	-	.	gene_id \"ENSG00000142606.16\"; gene_type \"protein_coding\"; gene_name \"MMEL1\"; level 2; hgnc_id \"HGNC:14668\"; havana_gene \"OTTHUMG00000000846.5\";\n");
        tao::pegtl::memory_input in(l, "string");
        scindo::gtf::detail::state s;
        REQUIRE( tao::pegtl::parse<scindo::gtf::detail::row, scindo::gtf::detail::action>(in, s) );
    }
}
TEST_CASE("parse tiny gtf", "[gtf]") {
    scindo::gtf::gtf_file G("data/tiny.gtf");

    G.parse([&](const std::string& p_seqname, const std::string& p_source, const std::string& p_feature, 
                const int64_t& p_start, const int64_t& p_end,
                const std::optional<double>& p_score, const char& p_strand, const std::optional<int>& p_frame,
                const std::vector<scindo::gtf::attribute>& p_attrs) {
    });
}
