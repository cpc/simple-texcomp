/* Common definitions for BCn algorithms
 */

#ifndef SIMPLE_MATHLIB_HPP
#define SIMPLE_MATHLIB_HPP

#include<cstdint>
#include<cstdio>

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

    // inline Vec3u8 operator-(uint8_t a) const
    // {
    //     return Vec3u8 {
    //         static_cast<uint8_t>(x - a),
    //         static_cast<uint8_t>(y - a),
    //         static_cast<uint8_t>(z - a),
    //     };
    // }
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

}  // namespace simple

#endif // SIMPLE_MATHLIB_HPP
