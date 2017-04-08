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
void bench(const char* func_name, calc_check_digit_f f, const std::vector<std::string>& inputs) {
	namespace ch = std::chrono;
	using namespace std::string_literals;
	using hc = ch::high_resolution_clock;
	const auto t0 = hc::now();
	[[gnu::unused]] std::uint8_t dst;
	for (auto&& i : inputs) dst = f(i);
	const auto t1 = hc::now();
	const auto t = t1 - t0;
	std::cout
		<< func_name << " : test::"
		<< ch::duration_cast<ch::milliseconds>(t).count() << "[ms] ("
		<< ch::duration_cast<ch::nanoseconds>(t).count() << "[ns])" << std::endl;
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


//@MaverickTse
//https://gist.github.com/MaverickTse/b78eff8fcc70962e0ee7a21b985bbaa9

/// Function to calculate check digit with AVX2 intrinsics
/// parm: a 11 digit number as string
/// ret: an unsigned integer 0~11
std::uint8_t get_check_digit_avx2(const std::string& query)
{
	if (query.length() != 11)
	{
		std::cerr << "Query string is not 11 characters" << std::endl;
		return 0xFF;
	}
	unsigned long long as_value{ 0 };
	try
	{
		as_value = std::stoull(query);
	}
	catch (std::invalid_argument e)
	{
		std::cerr << "Query string is not a valid number" << std::endl;
		return 0xFF;
	}
	catch (std::exception all)
	{
		std::cerr << "Error converting string to number." << std::endl;
		std::cerr << all.what() << std::endl;
		return 0xFF;
	}
	std::array<std::uint8_t, 32> digits{ 0 };

	std::array<short, 16> simd_result{ 0 }; // the 16bit intermediate results from SIMD

	for (int i = 0; i < 11; ++i)
	{
		digits[i] = as_value % 10;
		as_value /= 10;
	}
	// Note: the least significant digit is now the 1st item of array

	// Fire up SIMD for calculation

	// load P from array
	__m256i vP = _mm256_loadu_si256(reinterpret_cast<const __m256i*> (digits.data()));
	// Set Q, beware of param order
	__m256i vQ = _mm256_set_epi8(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 5, 4, 3, 2, 7, 6, 5, 4, 3, 2);
	// Multiply-add vP and vQ
	__m256i vR = _mm256_maddubs_epi16(vP, vQ);
	// Store vR
	_mm256_storeu_si256(reinterpret_cast<__m256i*>(simd_result.data()), vR);
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

/// Function to calculate check digit with SSE3 intrinsics
/// parm: a 11 digit number as string
/// ret: an unsigned integer 0~11
std::uint8_t get_check_digit_ssse3(const std::string& query)
{
	if (query.length() != 11)
	{
		std::cerr << "Query string is not 11 characters" << std::endl;
		return 0xFF;
	}
	unsigned long long as_value{ 0 };
	try
	{
		as_value = std::stoull(query);
	}
	catch (std::invalid_argument e)
	{
		std::cerr << "Query string is not a valid number" << std::endl;
		return 0xFF;
	}
	catch (std::exception all)
	{
		std::cerr << "Error converting string to number." << std::endl;
		std::cerr << all.what() << std::endl;
		return 0xFF;
	}
	std::array<std::uint8_t, 16> digits{ 0 };
	std::array<short, 16> simd_result{ 0 }; // the 16bit intermediate results from SIMD

	for (int i = 0; i < 11; ++i)
	{
		digits[i] = as_value % 10;
		as_value /= 10;
	}
	// Note: the least significant digit is now the 1st item of array

	// Fire up SIMD for calculation

	// load P from array
	__m128i vP = _mm_loadu_si128(reinterpret_cast<const __m128i*> (digits.data()));

	// Set Q, beware of order
	__m128i vQ = _mm_set_epi8(0, 0, 0, 0, 0, 6, 5, 4, 3, 2, 7, 6, 5, 4, 3, 2);
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

		bench("calc_check_digit_yumetodo_kai", calc_check_digit_yumetodo_kai, inputs);
		bench("calc_check_digit_ryogaelbtn", calc_check_digit_ryogaelbtn, inputs);
		bench("calc_check_digit_yumetodo", calc_check_digit_yumetodo, inputs);
		bench("calc_check_digit_ryogaelbtn2", calc_check_digit_ryogaelbtn2, inputs);
		bench("calc_check_digit_yumetodo_original", calc_check_digit_yumetodo_original, inputs);
		bench("get_check_digit_avx2", get_check_digit_avx2, inputs);
		bench("get_check_digit_ssse3", get_check_digit_ssse3, inputs);
		std::cout << "benchmark finish!" << std::endl;
	}
	catch (const std::exception& er) {
		std::cerr << er.what() << std::endl;
	}
}