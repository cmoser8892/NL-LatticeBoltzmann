#include <immintrin.h>
#include <gtest/gtest.h>

TEST(Extra_tests,avx) {
    __m256 A{};
    __m256 B{};
    __m256 AB = _mm256_add_ps(A,B);
    for(int i = 0; i < 8; i++){
        EXPECT_EQ(AB[i],0);
    }
}