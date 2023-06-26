#include "image_converter.h"

#include <immintrin.h> // intel instructions
#include <gtest/gtest.h>

/**
 * Checks if vectorization works and is enabled.
 * @test
 */
TEST(Extra_tests,avx) {
    __m256 A{};
    __m256 B{};
    __m256 AB = _mm256_add_ps(A,B);
    for(int i = 0; i < 8; i++){
        EXPECT_EQ(AB[i],0);
    }
}

/**
 * Checks weather or not the bmp headers are packed correctly and have the right sizes in memory.
 * @test
 */
TEST(Extra_tests, packed_and_counted_bmp_headers) {
    BMPFileHeader_t t;
    EXPECT_EQ(sizeof(t),14);
    BMPInfoHeader_t k;
    EXPECT_EQ(sizeof(k),40);
    BMPColorTable_t s;
    EXPECT_EQ(sizeof(s),20+16*4);
}