#include <iostream>
#include <cstdint>

int main() {
    // differences
    uint16_t a = 0x7536;
    uint16_t b = 0x8300;
    uint16_t c = 0x1812;

    // mask
    uint16_t u = 0x0030;
    uint16_t v = 0x1800;
    uint16_t w = 0x0000;

    // guess bit sum
    // Expected parity of the linear mask expression.
    uint16_t s = 0;

    // check
    uint64_t count = 0; // Number of assignments that violate the expected parity.
    // Wider counters include 0xffff without wrapping at the loop boundary.
    for (uint32_t x_value = 0; x_value <= UINT16_MAX; ++x_value) {
        for (uint32_t y_value = 0; y_value <= UINT16_MAX; ++y_value) {
            const uint16_t x = static_cast<uint16_t>(x_value);
            const uint16_t y = static_cast<uint16_t>(y_value);
            const uint16_t z = static_cast<uint16_t>(x + y);

            if ((uint16_t)(z ^ c) == (uint16_t)((x ^ a) + (y ^ b))) {
                uint16_t tmp = __builtin_popcount(x & u)
                             + __builtin_popcount(y & v)
                             + __builtin_popcount(z & w);
                if ((tmp & 1) != s) {
                    std::cout << (tmp & 1) << "\n";
                    count++;
                }
            }
        }
    }
    std::cout << "Total solutions: " << count << "\n";
    
    return 0;
}
