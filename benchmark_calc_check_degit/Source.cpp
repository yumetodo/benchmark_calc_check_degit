#include <iostream>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <array>
#include <random>
#include <chrono>
#include <random>
#include <vector>
#include <string>
#include <cstdint>
#include <algorithm>
#include <stdexcept>
#include <unordered_set>
#include <functional>
#include <limits>
#include <cctype>
#include <emmintrin.h>
#include <immintrin.h>
#define SPROUT_CONFIG_DISABLE_VARIABLE_TEMPLATES 1
#define SPROUT_CONFIG_FORCE_CXX14_CONSTEXPR 1
#include <sprout/array.hpp>
#if defined(_MSC_VER) && !defined(__c2__)
#pragma warning(disable: 5030)//警告	C5030	属性 'gnu::unused' は認識されません
#endif
std::mt19937 create_rand_engine() {
	std::random_device rnd;
	std::vector<std::uint_least32_t> v(10);// 初期化用ベクタ
	std::generate(v.begin(), v.end(), std::ref(rnd));// ベクタの初期化
	std::seed_seq seed(v.begin(), v.end());
	return std::mt19937(seed);// 乱数エンジン
}
std::string generate_input() {
	static auto mt = create_rand_engine();
	static std::uniform_int_distribution<std::uint16_t> dist(0u, 9u);
	std::string re;
	re.resize(11);
	std::generate(re.begin(), re.end(), []() -> char {
		return char('0' + dist(mt));
	});
	return re;
}
using calc_check_digit_f = std::uint8_t(*)(const std::string&);
struct test_result {
	std::string test_case;
	std::uint8_t expected;
	std::uint8_t actual;
};
std::ostream& operator<<(std::ostream& os, const test_result& t) {
	os << "testcase:" << t.test_case << ": expected(" << std::uint16_t(t.expected) << ") actual(" << std::uint16_t(t.actual) << ')';
	return os;
}
auto test(calc_check_digit_f f) {
	static const std::pair<std::string, std::uint8_t> testcaes[] = {
		{ "12345678901", 8 },
		{ "56661137362", 0 },
		{ "61671451309", 6 },
		{ "66383747390", 9 },
		{ "08065140466", 9 },
		{ "15473503678", 2 },
		{ "40113376378", 8 },
		{ "12226480680", 0 },
		{ "82046894873", 4 },
		{ "48880517043", 6 },
		{ "97816786786", 3 }
	};
	std::vector<test_result> fail;
	fail.reserve(sizeof(testcaes) / sizeof(*testcaes));
	for (auto&& t : testcaes) {
		try {
			const auto re = f(t.first);
			if (re != t.second) fail.push_back({ t.first, t.second, re });
		}
		catch (...) {
			fail.push_back({ t.first, t.second, std::uint8_t(0xff) });
		}
	}
	return fail;
}
void bench(const char* func_name, calc_check_digit_f f, const std::vector<std::string>& inputs) {
	namespace ch = std::chrono;
	using namespace std::string_literals;
	using hc = ch::high_resolution_clock;
	const auto test_re = test(f);
	auto v2s = [](const std::vector<test_result>& a) {
		std::stringstream ss;
		for (auto&& t : a) ss << t << std::endl;
		return ss.str();
	};
	const auto t0 = hc::now();
	[[gnu::unused]] std::uint8_t dst;
	for (auto&& i : inputs) dst = f(i);
	const auto t1 = hc::now();
	const auto t = t1 - t0;
	std::cout
		<< func_name << " : test::" << (test_re.empty() ? "pass:" : "fail:")
		<< ch::duration_cast<ch::milliseconds>(t).count() << "[ms] ("
		<< ch::duration_cast<ch::nanoseconds>(t).count() << "[ns])" << std::endl;
	if(!test_re.empty()) std::cout << v2s(test_re) << std::endl;
}
constexpr auto make_qn() {
	alignas(16) sprout::array<std::uint16_t, sizeof(__m256i) / sizeof(std::uint16_t)> re{};
	for (std::uint8_t i = 0, n = 1; i < re.size(); ++i, ++n) re[i] = (n < 7) ? n + 1 : n - 5;
	return re;
}
constexpr sprout::array<std::uint8_t, 1000> make_mod_table_ysr() {
	sprout::array<std::uint8_t, 1000> re{};
	for (size_t i = 0; i < 1000; ++i) {
		re[i] = 11 - (i % 11);
	}
	return re;
}
constexpr sprout::array<std::uint8_t, 1000> make_mod_table_yumetodo() {
	sprout::array<std::uint8_t, 1000> re{};
	for (size_t i = 0; i < 1000; ++i) {
		re[i] = i % 11;
	}
	return re;
}

std::uint8_t calc_check_digit_yumetodo_kai_simd(const std::string& n) noexcept(false) {
	static constexpr auto mod_table = make_mod_table_yumetodo();
	constexpr std::size_t num_of_digits = 11;
	if (num_of_digits != n.size()) throw std::runtime_error("n.digit must be 11");
	for (auto e : n) if (e < '0' || '9' < e) { throw std::runtime_error("in function calc_check_digit_yumetodo : iregal charactor detect.(" + n + ')'); }
	alignas(16) constexpr auto qn = make_qn();//0-7
	alignas(16) std::uint16_t n1[sizeof(__m256i) / sizeof(std::uint16_t)];
	for (std::size_t i = 0; i < num_of_digits; ++i) n1[i] = std::uint16_t(n[num_of_digits - 1 - i]);//reverse
	const __m256i pn1 = _mm256_sub_epi16(_mm256_load_si256(reinterpret_cast<const __m256i*>(n1)), _mm256_set1_epi16('0'));
	alignas(16) std::uint16_t tmp[sizeof(__m256i) / sizeof(std::uint16_t)];//0-63
	const __m256i qn1 = _mm256_load_si256(reinterpret_cast<const __m256i*>(qn.data()));
	const auto re = _mm256_mullo_epi16(pn1, qn1);
	_mm256_storeu_si256(reinterpret_cast<__m256i*>(tmp), re);
	std::uint16_t r = 0;
	for (std::size_t i = 0; i < num_of_digits; ++i) r += tmp[i];
	r = mod_table[r];
	return (0 == r || 1 == r) ? 0 : 11 - r;
}
std::uint8_t calc_check_digit_yumetodo_kai(const std::string& n) noexcept(false) {
	constexpr std::size_t num_of_digits = 11;
	if (num_of_digits != n.size()) throw std::runtime_error("n.digit must be 11");
	for (auto e : n) if (e < '0' || '9' < e) { throw std::runtime_error("in function calc_check_digit_yumetodo : iregal charactor detect.(" + n + ')'); }
	constexpr auto qn = make_qn();//0-7
	std::array<std::uint16_t, 11> pn;//0-9
	for (std::size_t i = 0; i < num_of_digits; ++i) pn[i] = std::uint16_t(n[num_of_digits - 1 - i] - '0');
	std::array<std::uint16_t, 11> tmp;//0-63
	for (std::size_t i = 0; i < num_of_digits; ++i) tmp[i] = pn[i] * qn[i];
	std::uint16_t r = 0;
	for (std::size_t i = 0; i < num_of_digits; ++i) r += tmp[i];
	r %= 11;
	return (0 == r || 1 == r) ? 0 : 11 - r;
}
std::uint8_t calc_check_digit_yumetodo(const std::string& n) noexcept(false) {
	if (11 != n.size()) throw std::runtime_error("n.digit must be 11");
	for (auto e : n) if (e < '0' || '9' < e) { throw std::runtime_error("in function calc_check_digit_yumetodo : iregal charactor detect.(" + n + ')'); }
	const std::uint8_t r = std::accumulate(n.rbegin(), n.rend(), std::pair<int, int>{}, [](const auto& s, const char& e) -> std::pair<int, int> {
		return { s.first + (e - '0') * ((5 < s.second) ? s.second - 4 : s.second + 2), s.second + 1 };
	}).first % 11;
	return (0 == r || 1 == r) ? 0 : 11 - r;
}
std::uint8_t calc_check_digit_yumetodo_original(const std::string& n) noexcept(false) {
	if (11 != n.size()) throw std::runtime_error("n.digit must be 11");
	const int r = std::accumulate(n.rbegin(), n.rend(), std::pair<int, int>{}, [](const auto& s, const char& e) -> std::pair<int, int> {
		if (!std::isdigit(e)) throw std::runtime_error("n.digit must be 11");
		return { s.first + (e - '0') * ((5 < s.second) ? s.second - 4 : s.second + 2), s.second + 1 };
	}).first % 11;
	return (0 == r || 1 == r) ? 0 : 11 - r;
}

//@proelbtn

constexpr std::uint16_t Q(std::uint8_t n) {
	return (1 <= n && n <= 6) ? n + 1 : n - 5;
}
std::uint8_t calc_check_digit_ryogaelbtn(const std::string& P) {
	if (11 != P.size()) throw std::runtime_error("P.digit must be 11");
	for (auto e : P) if (e < '0' || '9' < e) { throw std::runtime_error("in function calc_check_digit_ryogaelbtn : iregal charactor detect.(" + P + ')'); }
	unsigned int sum = 0;
	sum += static_cast<std::uint16_t>(P[10]) * Q(1);
	sum += static_cast<std::uint16_t>(P[9]) * Q(2);
	sum += static_cast<std::uint16_t>(P[8]) * Q(3);
	sum += static_cast<std::uint16_t>(P[7]) * Q(4);
	sum += static_cast<std::uint16_t>(P[6]) * Q(5);
	sum += static_cast<std::uint16_t>(P[5]) * Q(6);
	sum += static_cast<std::uint16_t>(P[4]) * Q(7);
	sum += static_cast<std::uint16_t>(P[3]) * Q(8);
	sum += static_cast<std::uint16_t>(P[2]) * Q(9);
	sum += static_cast<std::uint16_t>(P[1]) * Q(10);
	sum += static_cast<std::uint16_t>(P[0]) * Q(11);

	return sum % 11 <= 1 ? 0 : 11 - sum % 11;
}
std::uint8_t calc_check_digit_ryogaelbtn2(const std::string& P) {
	if (11 != P.size()) throw std::runtime_error("P.digit must be 11");
	try { std::stoull(P); }
	catch (...) { throw std::runtime_error("in function calc_check_digit_ryogaelbtn : iregal charactor detect.(" + P + ')'); }
	unsigned int sum = 0;
	sum += static_cast<std::uint16_t>(P[10]) * Q(1);
	sum += static_cast<std::uint16_t>(P[9]) * Q(2);
	sum += static_cast<std::uint16_t>(P[8]) * Q(3);
	sum += static_cast<std::uint16_t>(P[7]) * Q(4);
	sum += static_cast<std::uint16_t>(P[6]) * Q(5);
	sum += static_cast<std::uint16_t>(P[5]) * Q(6);
	sum += static_cast<std::uint16_t>(P[4]) * Q(7);
	sum += static_cast<std::uint16_t>(P[3]) * Q(8);
	sum += static_cast<std::uint16_t>(P[2]) * Q(9);
	sum += static_cast<std::uint16_t>(P[1]) * Q(10);
	sum += static_cast<std::uint16_t>(P[0]) * Q(11);

	return sum % 11 <= 1 ? 0 : 11 - sum % 11;
}
__m128i mullo_epi8(__m128i a, __m128i b)
{
	// unpack and multiply
	__m128i dst_even = _mm_mullo_epi16(a, b);
	__m128i dst_odd = _mm_mullo_epi16(_mm_srli_epi16(a, 8), _mm_srli_epi16(b, 8));
	// repack
#ifdef __AVX2__
	// only faster if have access to VPBROADCASTW
	return _mm_or_si128(_mm_slli_epi16(dst_odd, 8), _mm_and_si128(dst_even, _mm_set1_epi16(0xFF)));
#else
	return _mm_or_si128(_mm_slli_epi16(dst_odd, 8), _mm_srli_epi16(_mm_slli_epi16(dst_even, 8), 8));
#endif
}

//@YSRKEN
//http://qiita.com/YSRKEN/items/4ca7229c98640a71bdad

std::uint8_t calc_check_digit_ysrken(const std::string& str) noexcept(false) {
	static constexpr auto mod_table = make_mod_table_ysr();
	// 自分で作った文字列に対して入力チェックが必要なのかしら……？
	if (11 != str.size()) throw std::runtime_error("str.digit must be 11");
	for (auto e : str) if (e < '0' || '9' < e) { throw std::runtime_error("in function calc_check_digit_ysrken : iregal charactor detect.(" + str + ')'); }
	// __m128i型にマッピングし、'0'でマイナスすることで整数化する
	// 必然的に後ろ8ビット×5個=40ビット分はゴミが入ることになる
	__m128i p_n = _mm_loadu_si128(reinterpret_cast<const __m128i*>(str.c_str()));
	static const __m128i sub = _mm_set_epi8('0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0');
	p_n = _mm_sub_epi8(p_n, sub);
	// q_nはどうせ定数なので決め打ちする
	// p_nは反転処理すらしてないので注意
	static const __m128i q_n = _mm_set_epi8(0, 0, 0, 0, 0, 2, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6);
	// p_nとq_nとの掛け算
	const __m128i mul_pq = mullo_epi8(p_n, q_n);
	// 総和を計算する
	__m128i y = _mm_sad_epu8(mul_pq, _mm_setzero_si128());
	y = _mm_add_epi16(y, _mm_srli_si128(y, 8));
	return mod_table[_mm_cvtsi128_si32(y)];
}

//@MaverickTse
//https://gist.github.com/MaverickTse/b78eff8fcc70962e0ee7a21b985bbaa9

/// Function to calculate check digit with SSSE3 intrinsics
/// parm: a 11 digit number as string
/// ret: an unsigned integer 0~11
std::uint8_t get_check_digit_ssse3(const std::string& query)
{

	unsigned long long as_value{ 0 };
	std::array<short, 16> simd_result{ 0 }; // the 16bit intermediate results from SIMD
	if (11 != query.length()) throw std::runtime_error("str.digit must be 11");
	for (auto e : query) if (e < '0' || '9' < e) { throw std::runtime_error("in function get_check_digit_ssse3 : illegal character detected.(" + query + ')'); }
	__m128i vP = _mm_loadu_si128(reinterpret_cast<const __m128i*> (query.c_str()));

	__m128i vzero = _mm_set1_epi8('0');
	
	vP = _mm_sub_epi8(vP, vzero);

	// Set Q, beware of order
	__m128i vQ = _mm_set_epi8(0, 0, 0, 0, 0, 2, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6);
		
	// load P from array
	__m128i vP = _mm_loadu_si128(reinterpret_cast<const __m128i*> (digits.data()));

	// Multiply-add vP and vQ
	__m128i vR = _mm_maddubs_epi16(vP, vQ);
	// Store vR
	_mm_storeu_si128(reinterpret_cast<__m128i*>(simd_result.data()), vR);
	// our result
	int result{ 0 };
	for (int i = 0; i < 6; ++i)
	{
		result += simd_result[i];
	}

	result %= 11;
	if (result <= 1)
	{
		result = 0;
	}
	result = 11 - result;
	return static_cast<std::uint8_t>(result);
}
int main() {
	constexpr int test_times = 14000000;
	try {
		std::cout << "generating inputs..." << std::flush;
		std::vector<std::string> inputs(test_times);
		std::generate(inputs.begin(), inputs.end(), []() { return generate_input(); });
		std::cout
			<< "done." << std::endl
			<< "start benchmark mark:" << std::endl;

		bench("calc_check_digit_yumetodo_kai_simd", calc_check_digit_yumetodo_kai_simd, inputs);
		bench("calc_check_digit_yumetodo_kai", calc_check_digit_yumetodo_kai, inputs);
		bench("calc_check_digit_ryogaelbtn", calc_check_digit_ryogaelbtn, inputs);
		bench("calc_check_digit_yumetodo", calc_check_digit_yumetodo, inputs);
		bench("calc_check_digit_ryogaelbtn2", calc_check_digit_ryogaelbtn2, inputs);
		bench("calc_check_digit_yumetodo_original", calc_check_digit_yumetodo_original, inputs);
		bench("calc_check_digit_ysrken", calc_check_digit_ysrken, inputs);
		bench("get_check_digit_ssse3", get_check_digit_ssse3, inputs);
		std::cout << "benchmark finish!" << std::endl;
	}
	catch (const std::exception& er) {
		std::cerr << er.what() << std::endl;
	}
}
