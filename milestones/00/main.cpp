/**
 * 00 is general testing ground, quick and dirty
 */


#include <cstdint>

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
    WAVEFORMAT w;
    (void)w;
}