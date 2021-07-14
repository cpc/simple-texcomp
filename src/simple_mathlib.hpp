/* Common definitions for BCn algorithms
 */

#ifndef SIMPLE_MATHLIB_HPP
#define SIMPLE_MATHLIB_HPP

#include <cstdint>
#include <cstdio>

#ifndef NDEBUG
#include <cstring>
#endif

#include "simple_texcomp.hpp"

// TODO: It's growing out of hand, should be templated

namespace simple {

/* Endpoint interpolation constants */
constexpr decimal EP_LERP1 = F(1.0) / F(3.0);
constexpr decimal EP_LERP2 = F(2.0) / F(3.0);

/* Helper structs and math */
struct Vec3f
{
    decimal x;
    decimal y;
    decimal z;

    inline Vec3f operator+(const Vec3f &other) const
    {
        return Vec3f {
            x + other.x,
            y + other.y,
            z + other.z,
        };
    }

    inline Vec3f operator+(decimal a) const
    {
        return Vec3f {
            x + a,
            y + a,
            z + a,
        };
    }

    inline Vec3f operator-(const Vec3f &other) const
    {
        return Vec3f {
            x - other.x,
            y - other.y,
            z - other.z,
        };
    }

    inline Vec3f operator-(decimal a) const
    {
        return Vec3f {
            x - a,
            y - a,
            z - a,
        };
    }

    inline Vec3f operator*(decimal a) const
    {
        return Vec3f {
            x * a,
            y * a,
            z * a,
        };
    }

    inline Vec3f operator/(decimal a) const
    {
        return Vec3f {
            x / a,
            y / a,
            z / a,
        };
    }

    inline decimal dot(const Vec3f &other) const
    {
        // return std::fmaf(x, other.x, std::fmaf(y, other.y, fmaf(z, other.z, F(0.0))));
        return (x * other.x) + (y * other.y) + (z * other.z);
    }
};

struct Vec2f
{
    decimal x;
    decimal y;

    inline Vec2f operator+(const Vec2f &other) const
    {
        return Vec2f {
            x + other.x,
            y + other.y,
        };
    }

    inline Vec2f operator+(decimal a) const
    {
        return Vec2f {
            x + a,
            y + a,
        };
    }

    inline Vec2f operator-(const Vec2f &other) const
    {
        return Vec2f {
            x - other.x,
            y - other.y,
        };
    }

    inline Vec2f operator-(decimal a) const
    {
        return Vec2f {
            x - a,
            y - a,
        };
    }

    inline Vec2f operator*(decimal a) const
    {
        return Vec2f {
            x * a,
            y * a,
        };
    }

    inline Vec2f operator/(decimal a) const
    {
        return Vec2f {
            x / a,
            y / a,
        };
    }

    inline decimal dot(const Vec2f &other) const
    {
        return (x * other.x) + (y * other.y);
    }
};

struct Vec3i
{
    int x;
    int y;
    int z;

    inline Vec3i operator>>(unsigned int nbits) const
    {
        return Vec3i {
            x >> nbits,
            y >> nbits,
            z >> nbits,
        };
    }

    inline bool operator==(const Vec3i &other) const
    {
        return (x == other.x) && (y == other.y) && (z == other.z);
    }
};

struct Vec3u32
{
    uint32_t x;
    uint32_t y;
    uint32_t z;
};

struct Vec3u8
{
    uint8_t x;
    uint8_t y;
    uint8_t z;

    inline Vec3u8 operator>>(uint8_t nbits) const
    {
        return Vec3u8 {
            static_cast<uint8_t>(x >> nbits),
            static_cast<uint8_t>(y >> nbits),
            static_cast<uint8_t>(z >> nbits),
        };
    }

    inline bool operator==(const Vec3u8 &other) const
    {
        return (x == other.x) && (y == other.y) && (z == other.z);
    }

    inline Vec3u8 operator-(const Vec3u8 &other) const
    {
        // can and will overflow
        return Vec3u8 {
            static_cast<uint8_t>(x - other.x),
            static_cast<uint8_t>(y - other.y),
            static_cast<uint8_t>(z - other.z),
        };
    }

    inline uint16_t dot(const Vec3u8 &other) const
    {
        // Arithmetic dot product, returns the 16 LSB of the result
        uint32_t xx = x * other.x;
        uint32_t yy = y * other.y;
        uint32_t zz = z * other.z;

        uint32_t res32 = xx + yy + zz;
        uint16_t res16 = (uint16_t)(res32);

        if (res32 != (uint32_t)(res16))
        {
            printf("WARNING: Dot product overflow! (%3d %3d %3d).(%3d %3d %3d)\n",
                x, y, z, other.x, other.y, other.z);
        }

        return (uint16_t)(xx + yy + zz);
    }

    inline uint16_t dot16(const Vec3u8 &other) const
    {
        // Dot product assuming fixed precision
        // Keeps the multiplication precision but rounds the bits added by +
        uint32_t xx = x * other.x;
        uint32_t yy = y * other.y;
        uint32_t zz = z * other.z;

        // max. possible value is 255*255*3 = 0b10.1111101000000011 (Q18.16)
        return (uint16_t)((xx + yy + zz) >> 2);  // round off the two additions
    }

    inline uint16_t sataccdot(const Vec3u8 &other, uint32_t acc) const
    {
        // Saturating dot product (overflow clamps to 0xffff)
        uint32_t xx = (uint32_t)x * (uint32_t)other.x;
        uint32_t yy = (uint32_t)y * (uint32_t)other.y;
        uint32_t zz = (uint32_t)z * (uint32_t)other.z;

        uint32_t res32 = acc + xx + yy + zz;

        if (res32 > 0xffff)
        {
            return 0xffff;
        }

        return (uint16_t)(res32);
    }

    inline uint32_t dot32(const Vec3u8 &other) const
    {
        // Dot product assuming fixed precision
        // Keeps the multiplication precision but rounds the bits added by +
        uint32_t xx = x * other.x;
        uint32_t yy = y * other.y;
        uint32_t zz = z * other.z;

        // max. possible value is 255*255*3 = 0b10.1111101000000011 (Q18.16)
        return (uint32_t)(xx + yy + zz);  // return the full result
    }
};

struct Vec3u16
{
    uint16_t x;
    uint16_t y;
    uint16_t z;

    inline uint16_t dotfi(const Vec3u16 &other) const
    {
        // Dot product assuming fixed precision
        // Keeps the multiplication precision but rounds the bits added by +
        uint32_t xx = x * other.x;
        uint32_t yy = y * other.y;
        uint32_t zz = z * other.z;

        return (uint16_t)((xx + yy + zz) >> 2);  // round off the two additions
    }
};

inline decimal fclamp(decimal a, decimal amin, decimal amax)
{
    const decimal min = a < amin ? amin : a;
    return min > amax ? amax : min;
}

inline int iclamp(int a, int amin, int amax)
{
    const int min = a < amin ? amin : a;
    return min > amax ? amax : min;
}

inline uint8_t u8clamp(uint8_t a, uint8_t amin, uint8_t amax)
{
    const uint8_t min = a < amin ? amin : a;
    return min > amax ? amax : min;
}

inline decimal fmin(decimal a, decimal b)
{
    if (a < b)
    {
        return a;
    } else
    {
        return b;
    }
}

inline decimal fmax(decimal a, decimal b)
{
    if (a > b)
    {
        return a;
    } else
    {
        return b;
    }
}

inline uint8_t u8min(uint8_t a, uint8_t b)
{
    if (a < b)
    {
        return a;
    } else
    {
        return b;
    }
}

inline uint8_t u8max(uint8_t a, uint8_t b)
{
    if (a > b)
    {
        return a;
    } else
    {
        return b;
    }
}

inline Vec3f min3f(const Vec3f &a, const Vec3f &b)
{
    return Vec3f {
        (decimal)fmin(a.x, b.x),
        (decimal)fmin(a.y, b.y),
        (decimal)fmin(a.z, b.z),
    };
}

inline Vec3f max3f(const Vec3f &a, const Vec3f &b)
{
    return Vec3f {
        (decimal)fmax(a.x, b.x),
        (decimal)fmax(a.y, b.y),
        (decimal)fmax(a.z, b.z),
    };
}

inline Vec3f clamp3f(const Vec3f &a, decimal amin, decimal amax)
{
    return Vec3f {
        fclamp(a.x, amin, amax),
        fclamp(a.y, amin, amax),
        fclamp(a.z, amin, amax),
    };
}

inline Vec3i clamp3i(const Vec3i &a, int amin, int amax)
{
    return Vec3i {
        iclamp(a.x, amin, amax),
        iclamp(a.y, amin, amax),
        iclamp(a.z, amin, amax),
    };
}

inline Vec3u8 clamp3u8(const Vec3u8 &a, uint8_t amin, uint8_t amax)
{
    return Vec3u8 {
        u8clamp(a.x, amin, amax),
        u8clamp(a.y, amin, amax),
        u8clamp(a.z, amin, amax),
    };
}

inline Vec2f clamp2f(const Vec2f &a, decimal amin, decimal amax)
{
    return Vec2f {
        fclamp(a.x, amin, amax),
        fclamp(a.y, amin, amax),
    };
}

inline decimal fabs(decimal a)
{
    if (a < 0) {
        return -a;
    }
    else
    {
        return a;
    }
}

inline Vec2f abs2f(const Vec2f &a)
{
    return Vec2f { (decimal)fabs(a.x), (decimal)fabs(a.y) };
}

/* Compute squared distance of two vectors */
inline decimal distsq3f(const Vec3f &a, const Vec3f &b)
{
    Vec3f diff = a - b;
    return diff.dot(diff);
}

inline decimal distsq2f(const Vec2f &a, const Vec2f &b)
{
    Vec2f diff = a - b;
    return diff.dot(diff);
}

/* Round floating point number and return it as unsigned integer */
inline uint32_t froundu(decimal a)
{
    return (uint32_t)(a + F(0.5));
}

inline Vec3u8 satadd(const Vec3u8 &inp, const Vec3u8 &other)
{
    uint8_t x = inp.x + other.x;
    uint8_t y = inp.y + other.y;
    uint8_t z = inp.z + other.z;

    if (x < inp.x)
    {
        x = 255;
    }

    if (y < inp.y)
    {
        y = 255;
    }

    if (z < inp.z)
    {
        z = 255;
    }

    return Vec3u8 { x, y, z };
}

inline Vec3u8 satsub(const Vec3u8 &inp, const Vec3u8 &other)
{
    uint8_t x = inp.x - other.x;
    uint8_t y = inp.y - other.y;
    uint8_t z = inp.z - other.z;

    if (x > inp.x)
    {
        x = 0;
    }

    if (y > inp.y)
    {
        y = 0;
    }

    if (z > inp.z)
    {
        z = 0;
    }

    return Vec3u8 { x, y, z };
}

template<typename T, typename U>
inline T rshift_round(T x, U n)
{
    return (x >> n) + ((x >> (n - 1)) & 1);
}

/* Debug prints */
#ifndef NDEBUG

inline void print_bin(unsigned int num, unsigned int nb)
{
	for (unsigned int b = 0; b < nb; ++b)
	{
		unsigned int x = (num >> (nb-1-b)) & 1;
		printf("%d", x);
	}
}

inline void print_bin_(unsigned int num, unsigned int nb, unsigned int frac)
{
    unsigned int dec_point = nb - frac;
	for (unsigned int b = 0; b < nb; ++b)
	{
        if (b == dec_point)
        {
            printf(".");
        }
		unsigned int x = (num >> (nb-1-b)) & 1;
		printf("%d", x);
	}
}

inline void print_minmax(const char *pre, Vec3f mincol, Vec3f maxcol)
{
    Vec3f mincol2 = mincol * F(255.0);
    Vec3f maxcol2 = maxcol * F(255.0);

    printf("%s mincol: %5.3f %5.3f %5.3f  %3.0f %3.0f %3.0f\n",
        pre,
        (double)mincol.x, (double)mincol.y, (double)mincol.z,
        (double)mincol2.x, (double)mincol2.y, (double)mincol2.z
    );
    printf("    maxcol: %5.3f %5.3f %5.3f  %3.0f %3.0f %3.0f\n",
        (double)maxcol.x, (double)maxcol.y, (double)maxcol.z,
        (double)maxcol2.x, (double)maxcol2.y, (double)maxcol2.z
    );
}

inline void print_minmax(const char *pre, Vec3u8 mincol, Vec3u8 maxcol)
{
    int l = strlen(pre);
    printf("%*s mincol: %3d %3d %3d\n", l, pre, mincol.x, mincol.y, mincol.z);
    printf("%*s maxcol: %3d %3d %3d\n", l, "", maxcol.x, maxcol.y, maxcol.z);
}

inline void print_minmax(const char *pre, Vec3u8 mincol, Vec3u16 maxcol)
{
    int l = strlen(pre);
    printf("%*s mincol: %3d %3d %3d\n", l, pre, mincol.x, mincol.y, mincol.z);
    printf("%*s maxcol: %3d %3d %3d\n", l, "", maxcol.x, maxcol.y, maxcol.z);
}

inline double fixed_to_double(unsigned int x, unsigned int frac)
{
    return (double)(x) / ( (double)(1 << frac) );
}

inline void print_fixed(const char *pre, unsigned int x, unsigned int frac)
{
    int l = strlen(pre);
    printf("%*s %10u  ", l, pre, x);
    print_bin_(x, 32, frac);
    printf("  %13.8f\n", fixed_to_double(x, frac));
}

inline void print_fixed8(const char *pre, unsigned int x, unsigned int frac)
{
    int l = strlen(pre);
    double xfi = fixed_to_double(x, frac);

    printf("%*s %3u  ", l, pre, x);
    print_bin_(x, 8, frac);
    printf("  %12.8f", xfi);
}

#endif // NDEBUG

 /* Approximate 1/x in fixed point arithmetic.
  *
  * Assumes x is in Q2.16 format (unsigned). Returns Q10.22.
  */
inline uint32_t approx_inv32(uint32_t x)
{
    uint32_t xx = x;

    // First, scale the input to be within [0.5, 1.0]
    int32_t scale = (xx > 0xffff) << 4;
    xx >>= scale;

    uint32_t shift = (xx > 0xff) << 3;
    xx >>= shift;
    scale |= shift;

    shift = (xx > 0xf ) << 2;
    xx >>= shift;
    scale |= shift;

    shift = (xx > 0x3 ) << 1;
    xx >>= shift;
    scale |= shift;

    scale |= (xx >> 1);  // now, scale is log2(x)
    scale = 15 - scale;

    const uint32_t shl = (scale < 0) ?      0 : scale;
    const uint32_t shr = (scale < 0) ? -scale :     0;
    const uint32_t x_sc = (x << shl) >> (shr+1); // x_sc is 0.5--1.0, so Q0.16 -> Q0.15

    // Then, compute the initial estimate
    constexpr uint32_t A = 0x0000F0F1; // 32.0 / 17.0  Q2.15
    constexpr uint32_t B = 0xB4B4B4B5; // 48.0 / 17.0  Q2.30

    const uint32_t A_x_sc = A * x_sc;  // Q4.30 but x_sc is <= 1 so the two MSB are unused => Q2.30
    const uint32_t init = B - A_x_sc;  // Q2.30 - Q2.30 = Q3.30 -> Q2.30 (won't overflow)

    constexpr uint32_t ONE = (1 << 15) - 1;  // 1 in Q0.15

    // Newthon-Raphson iterations
    // 1st
    uint32_t y0 = init >> 15;             // Q2.15
    uint32_t y00 = init;                  // Q2.30
    uint32_t tmp = (x_sc * y0) >> 15;     // Q0.15 * Q2.15 = Q2.30 -> Q2.15
    uint32_t y1 = y00 + y0 * (ONE - tmp); // Q2.30 + (Q2.15 * (Q0.15 - Q2.15))

    // 2nd
    y0 = y1 >> 15;
    y00 = y1;
    tmp = (x_sc * y0) >> 15;
    y1 = y00 + y0 * (ONE - tmp);

    // 3rd
    y0 = y1 >> 15;
    y00 = y1;// >> 1;
    tmp = (x_sc * y0) >> 15;
    y1 = y00 + y0 * (ONE - tmp);

    // The result is scaled down now, we need to scale it back
    y1 >>= 8; // Q10.22
    return (y1 << shl) >> shr;
}

/* Compute log2 (i.e., the position of the highest bit set) */
inline uint32_t log2(uint32_t x)
{
    uint32_t scale = (x > 0xffff) << 4;
    x >>= scale;

    uint32_t shift = (x > 0xff) << 3;
    x >>= shift;
    scale |= shift;

    shift = (x > 0xf ) << 2;
    x >>= shift;
    scale |= shift;

    shift = (x > 0x3 ) << 1;
    x >>= shift;
    scale |= shift;

    scale |= (x >> 1);

    return scale;
}

}  // namespace simple

#endif // SIMPLE_MATHLIB_HPP
