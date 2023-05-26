/**
 * 00 is general testing ground, quick and dirty
 * when i dont want it in tests
 */


#include <cstdint>


#include <cassert>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>

typedef struct __attribute__((packed))
{
    uint16_t wFormatTag;
    uint16_t nChannels;
    uint32_t nSamplesPerSec;
    uint32_t nAvgBytesPerSec;
    uint16_t nBlockAlign;
}
WAVEFORMAT;

int main()
{
    // todo find out about unique ptrs
    WAVEFORMAT w;
    (void)w;
}