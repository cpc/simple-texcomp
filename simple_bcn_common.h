/* Common definitions for BCn algorithms
 */

#ifndef SIMPLE_BCN_COMMON_H
#define SIMPLE_BCN_COMMON_H

#include<stdint.h>
#include<stdio.h>

#include<cmath>

/* Number of channels used per pixel (each is assumed to be 8 bits) */
#define NCH_RGB  3

/* Endpoint interpolation constants */
#define EP_LERP1  1.0f / 3.0f
#define EP_LERP2  2.0f / 3.0f

/* Helper structs and math */
struct Vec3f
{
    double x;
    double y;
    double z;

    inline Vec3f operator+(const Vec3f &other) const
    {
        return Vec3f {
            x + other.x,
            y + other.y,
            z + other.z,
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

    inline Vec3f operator-(double a) const
    {
        return Vec3f {
            x - a,
            y - a,
            z - a,
        };
    }

    inline Vec3f operator*(double a) const
    {
        return Vec3f {
            x * a,
            y * a,
            z * a,
        };
    }

    inline Vec3f operator/(double a) const
    {
        return Vec3f {
            x / a,
            y / a,
            z / a,
        };
    }

    inline double dot(const Vec3f &other) const
    {
        // return std::fmaf(x, other.x, std::fmaf(y, other.y, fmaf(z, other.z, 0.0)));
        return (x * other.x) + (y * other.y) + (z * other.z);
    }
};

struct Vec2f
{
    double x;
    double y;

    inline Vec2f operator+(const Vec2f &other) const
    {
        return Vec2f {
            x + other.x,
            y + other.y,
        };
    }

    inline Vec2f operator+(double a) const
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

    inline Vec2f operator-(double a) const
    {
        return Vec2f {
            x - a,
            y - a,
        };
    }

    inline Vec2f operator*(double a) const
    {
        return Vec2f {
            x * a,
            y * a,
        };
    }

    inline Vec2f operator/(double a) const
    {
        return Vec2f {
            x / a,
            y / a,
        };
    }

    inline double dot(const Vec2f &other) const
    {
        return (x * other.x) + (y * other.y);
    }
};

// inline double minf(double a, double b)
// {
//     return (a < b) ? a : b;
// }

// inline double maxf(double a, double b)
// {
//     return (a > b) ? a : b;
// }

// inline double absf(double a)
// {
//     return (a < 0) ? -a : a;
// }

inline double fclamp(double a, double amin, double amax)
{
    const double min = a < amin ? amin : a;
    return min > amax ? amax : min;
}

inline Vec3f min3f(const Vec3f &a, const Vec3f &b)
{
    return Vec3f {
        fmin(a.x, b.x),
        fmin(a.y, b.y),
        fmin(a.z, b.z),
    };
}

inline Vec3f max3f(const Vec3f &a, const Vec3f &b)
{
    return Vec3f {
        fmax(a.x, b.x),
        fmax(a.y, b.y),
        fmax(a.z, b.z),
    };
}

inline Vec3f clamp3f(const Vec3f &a, double amin, double amax)
{
    return Vec3f {
        fclamp(a.x, amin, amax),
        fclamp(a.y, amin, amax),
        fclamp(a.z, amin, amax),
    };
}

inline Vec2f clamp2f(const Vec2f &a, double amin, double amax)
{
    return Vec2f {
        fclamp(a.x, amin, amax),
        fclamp(a.y, amin, amax),
    };
}

inline Vec2f abs2f(const Vec2f &a)
{
    return Vec2f { fabs(a.x), fabs(a.y) };
}

/* Compute squared distance of two vectors */
inline double distsq3f(const Vec3f &a, const Vec3f &b)
{
    Vec3f diff = a - b;
    return diff.dot(diff);
}

inline double distsq2f(const Vec2f &a, const Vec2f &b)
{
    Vec2f diff = a - b;
    return diff.dot(diff);
}

/* Find min/max color as a corners of a bounding box of the block */
inline void find_minmaxcolor_bbox(
    const Vec3f block[16],
    Vec3f *mincol,
    Vec3f *maxcol
){
    *mincol = { 1.0f, 1.0f, 1.0f };
    *maxcol = { 0.0f, 0.0f, 0.0f };

    for (int i = 0; i < 16; ++i)
    {
        *mincol = min3f(*mincol, block[i]);
        *maxcol = max3f(*maxcol, block[i]);
    }
}

#endif // SIMPLE_BCN_COMMON_H
