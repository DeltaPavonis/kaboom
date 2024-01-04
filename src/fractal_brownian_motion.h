#ifndef FRACTAL_BROWNIAN_MOTION_H
#define FRACTAL_BROWNIAN_MOTION_H

#include <numbers>
#include <cmath>
#include <numeric>
#include <array>
#include <random>
#include <optional>
#include "vec3d.h"

/* `seed` is the random seed used in generating the explosion. Each seed corresponds to exactly
one explosion. The user can choose to provide it or not, so I use `std::optional` for its type. */
std::optional<double> seed;
/* `GRID_CELL_SIZE` = The side length of all the cubes in the 3D grid we define for computing
Value Noise. */
constexpr double GRID_CELL_SIZE = 1.;
/* `CELL_HASH_VECTOR` is used to hash cubes in the 3D grid to a (theoretically) unique number;
the hash of a cube is equal to the dot product of its minimum-coordinate corner point and
`CELL_HASH_VECTOR`. We will need to hash cubes in `get_value`, which is the function responsible
for returning the values assigned to the intersection points in the 3D grid we define for
computing Value Noise. */
constexpr Vec3D CELL_HASH_VECTOR{
    /* This is the same vector that `ssloy` uses. I discuss choices of `CELL_HASH_VECTOR`
    in the comments for `get_value`.  */
    1, 57, 113
};

/* Sets `seed` to a random `double` in the range `[0, 1000000]`, if it does not have a value. */
void generate_seed_if_not_set() {
    if (!seed) {
        std::mt19937 generator{std::random_device{}()};
        std::uniform_real_distribution<> dist(0, 1000000);
        seed = dist(generator);
    }
}

/* Hashes the `double` `d` to a `double` in the range `[0, 1)`. */
auto get_hash(double d) {
    /* We use the same hash function in the tutorial; we return the fractional part of
    `sin(d) * seed` (this is where `seed` comes in). The fractional part by definition
    is in the range [0, 1), as desired. */
    auto temp = *seed * std::sin(d);
    return temp - std::floor(temp);
}

/* Returns the value at the point `p`, where `p` is an intersection point in the 3D grid we
use to compute the Value Noise. That is, `p` is expected to be the corner of some cube
in the 3D grid. This means the components of `p` should all be multiples of `GRID_CELL_SIZE`.

Each grid intersection point is assigned a value in the range [0, 1] (or well, more specifically,
[0, 1)). */
auto get_value(const Point3D &p) {
    /* As in the tutorial, the value of a grid intersection point p = (x, y, z) is given by
    the hash value of the dot product of p with our fixed `CELL_HASH_VECTOR`. This approach
    (using a hash function to generate random values) saves storage, because it allows us
    to forgo precomputing and storing the random values for every grid intersection point.

    The choice of `CELL_HASH_VECTOR` is important. We don't want `dot(p, CELL_HASH_VECTOR)`
    to equal the same thing for different points `p` (they shouldn't, because we're supposed
    to generate a random real number in [0, 1] for every grid intersection point, and random
    real numbers are theoretically distinct with probability 1). Thus, it would be wise to
    make the components of `CELL_HASH_VECTOR` three different irrational numbers, to minimize
    the probability of a collision between two different grid intersection points (I think).
    
    However, for now, I'm just using the same `CELL_HASH_VECTOR` as ssloy: {1, 57, 113}. The
    sole reason I'm doing this is because changing the `CELL_HASH_VECTOR` also changes the
    generated explosion (because changing `CELL_HASH_VECTOR` changes the values we give to
    each grid intersection point, which in turn affects the Value Noise and so changes the
    output of `fractal_brownian_motion`), and I haven't been able to find a random seed I
    like as much with a different `CELL_HASH_VECTOR`. */
    return get_hash(dot(p, CELL_HASH_VECTOR));
}

/* Computes and returns the 3D Value Noise (a real number in the range [0, 1]) at the point `p`. */
auto value_noise_at(const Point3D &p) {
    /* This function generates noise (real numbers in the range [0, 1]). ssloy's tutorial says we
    will "be helped by the Perlin's noise", but this is wrong. In actuality, we generate noise
    using the VALUE NOISE function. While Value Noise and Perlin Noise are closely related concepts
    (both generate noise given a single point `p`), they are conceptually and algorithmically
    different.
    
    N-dimensional PERLIN Noise works as follows. First, define an N-dimensional grid. Then, assign
    a fixed random unit gradient vector to every grid intersection point. Now, to find the Perlin
    Noise at a given point P = (x, y, z), do the following.
        1. Find the unique N-dimensional grid cell in which P lies. For example, in 3D Perlin Noise,
        we will need to find the unique cube in the grid in which P lies.
        2. Iterate through all 2^N vertices of that N-dimensional grid cell. For each vertex of
        that grid cell v, first compute (P - v); this is the displacement vector from that vertex
        to the given point P. Then, compute the dot product of that displacement vector with the
        given unit gradient vector for the vertex `v`. After this, we will be left with 2^N scalars,
        one for each vertex of the N-dimensional grid cell.
        3. Interpolate between the 2^N scalars in Part 2, based on the relative position of P within
        the grid cell that contains it. This will yield a single value, which after scaling it to
        make it fall in the range [0, 1], we will return.
            Note: A defining characteristic of Perlin Noise is that the noise returned from any
            grid intersection point will be 0. This is because if the given point is one of the
            grid intersection points, then its dot product with that grid intersection point will
            be 0, and furthermore, because it is located on that point, the interpolation between
            the corners of its square will bias completely towards the grid intersection point it
            is equal to, thus making the noise 0 at that point.
            Note: The range of Perlin Noise is not actually [0, 1]; it is [-sqrtN / 2, sqrtN / 2]
            where N is the number of dimensions used (a proof can be found at this link:
            (https://digitalfreepen.com/2017/06/20/range-perlin-noise.html). Note that this range
            also scales linearly with the side lengths of the grid cells; specifically, a grid
            cell side length of k results in the range [k * -sqrtN / 2, k * sqrtN / 2].

    N-dimensional VALUE Noise works as follows. First, define an N-dimensional grid. Then, assign
    a fixed random VALUE (not an unit gradient vector, unlike Perlin Noise) to every grid
    intersection point. Now, to find the Value Noise at a given point P = (x, y, z):
        1. Find the unique N-dimensional grid cell in which P lies. For example, in 3D Value Noise,
        we will need to find the unique cube in the grid in which P lies. This is identical to the
        first step in Perlin Noise.
        2. Consider all 2^N vertices of that N-dimensional grid cell (just as in Perlin Noise).
        3. The value noise at P is found by simply interpolating between the given values of those
        vertices, based on the relative position of P within the grid cell that contains it. This
        is different from Perlin Noise, where we interpolate between the dot products of the
        displacement vectors from each corner to P and the given gradient vectors at each
        corner.
        4. Finally, scale the interpolated value to your desired range.
            Note: In my code, I require that the fixed random value we assign to every grid
            intersection point at the beginning be in [0, 1]. As a result, the value we receive
            from step 3 must also be in the range [0, 1], because it is found by interpolating
            between the initial values given to the grid intersection points. Thus, in my
            implementation, we can skip this step entirely; we are guaranteed that the interpolated
            value falls in the desired range of [0, 1] already (because all the values of the grid
            intersection points fall in [0, 1]), so no scaling is needed.
    */

    /* Find the unique 3D grid cell (grid cube) in which `p` = (x, y, z) lies. Specifically,
    we will find the minimum-coordinate corner (the corner with the minimum x-, y-, and z-
    coordinates all at once, guaranteed to exist because the grid cells are parallel to the
    x-, y-, and z- axes) of the grid cube which contains `p`. Mathematically, this is equivalent
    to finding the unique grid intersection point (gx, gy, gz) where x/y/z are in the range
    [gx/gy/gz, gx/gy/gz + GRID_CELL_SIZE). We store the unique minimum-coordinate corner in the
    variable `min_coordinate_corner`. */
    Point3D min_coordinate_corner{
        /* Our 3D grid of cubes is defined as follows: a point (a, b, c) is a grid intersection
        point if and only if a, b, and c are all multiples of `GRID_CELL_SIZE`. Clearly, this
        defines a 3D grid of cubes with side length `GRID_CELL_SIZE`, which is exactly what we
        want.
        
        Then, the minimum-coordinate corner of the unique 3D grid cell in which `p` lies is simply
        the point where each coordinate is the largest multiple of `GRID_CELL_SIZE` that is at
        most the corresponding coordinate in `p`. Thus, the x/y/z coordinates of the
        minimum-coordinate corner of the cube are floor(x/y/z / GRID_CELL_SIZE) * GRID_CELL_SIZE. */
        std::floor(p.x / GRID_CELL_SIZE) * GRID_CELL_SIZE,
        std::floor(p.y / GRID_CELL_SIZE) * GRID_CELL_SIZE,
        std::floor(p.z / GRID_CELL_SIZE) * GRID_CELL_SIZE,
    };

    /* Now, the value noise at the point `p` is found by interpolating between the values of the
    grid cube's eight corners. First, for convenience, let's store all those values in
    `cell_corner_values`. */
    const std::array cell_corner_values{
        /* Bottom face of the cube (the four vertices with minimum z-coordinate). */
        get_value(min_coordinate_corner + GRID_CELL_SIZE * Point3D{0, 0, 0}),
        get_value(min_coordinate_corner + GRID_CELL_SIZE * Point3D{1, 0, 0}),
        get_value(min_coordinate_corner + GRID_CELL_SIZE * Point3D{0, 1, 0}),
        get_value(min_coordinate_corner + GRID_CELL_SIZE * Point3D{1, 1, 0}),
        /* Top face of the cube (the four vertices with maximum z-coordinate) */
        get_value(min_coordinate_corner + GRID_CELL_SIZE * Point3D{0, 0, 1}),
        get_value(min_coordinate_corner + GRID_CELL_SIZE * Point3D{1, 0, 1}),
        get_value(min_coordinate_corner + GRID_CELL_SIZE * Point3D{0, 1, 1}),
        get_value(min_coordinate_corner + GRID_CELL_SIZE * Point3D{1, 1, 1}),
    };
    /* ^^ ssloy's implementation works differently. They do calculate the hash of the point
    `min_coordinate_corner`, but then they hardcode the values passed to `get_hash`
    inside `get_value` for the other seven corners, using the observation that
    dot(min_coordinate_corner + Point3D{a, b, c}, CELL_HASH_VECTOR) is simply equal
    to dot(min_coordinate_corner, CELL_HASH_VECTOR) + a * CELL_HASH_VECTOR.x +
    b * CELL_HASH_VECTOR.y + c * CELL_HASH_VECTOR.z. */

    /* Now, compute the interpolation parameters, which we will use to interpolate between
    the values of the grid cube's eight corners. Here, we will just be using trilinear interpolation
    (repeated linear interpolation along the x-/y-/z-axes), so we need to find the corresponding
    interpolation parameters along the x-/y-/z- axes. To do this, simply observe that the
    interpolation parameter along a given axis is equivalent to the normalized value of `p`'s
    coordinate along that axis, with respect to the grid cube. For instance, if the cube's
    x-value ranges from A to B, then the interpolation parameter along the x-axis is given
    by (p.x - A) / (B - A). As a result, the interpolation parameters along the x/y/z-axes are
    determined by the RELATIVE POSITION of `p` inside the cube; the closer `p` is to a given
    corner, the greater the weight of that corner's value in determining the noise value at `p`.
    
    Finally, observe that the normalized position of `p` is equivalent to the offset of `p`
    from the minimum coordinate corner, divided by the side length of the cube (which is
    just `GRID_CELL_SIZE`). Thus, the interpolation parameters are succinctly computed as
    follows:*/
    auto interpolation_parameters = (p - min_coordinate_corner) / GRID_CELL_SIZE;
    {
        /* Now, we perform Hermitian Smoothing. This is a very common technique, where we
        apply a Hermite polynomial to each of the interpolation parameters. Commonly used
        polynomials include f(x) = 3x^2 - 2x^3 (called "smoothstep"). We use an "improved"
        version of Smoothstep, invented by Ken Perlin himself: 6x^5 - 15x^4 + 10x^3.
        
        The purpose of Hermitian Smoothing is to smooth the transition between points in
        the overall structure. From my testing, it seems to result in the generation of
        more structures on the surface of the explosion, giving it a better "cloud"-like
        look. */
        auto &[x, y, z] = interpolation_parameters;
        x = x * x * x * (x * (6 * x - 15) + 10);
        y = y * y * y * (y * (6 * y - 15) + 10);
        z = z * z * z * (z * (6 * z - 15) + 10);
    }

    /* Perform Trilinear Interpolation between the eight corners of the grid cube to find
    the Value Noise at the given point `p` (which is inside the grid cube), and return that.
    Trilinear interpolation is equivalent to three repeated linear interpolations, one along
    each of the axes. Here, I choose to first interpolate along the x-axis, then along the
    y-axis, and finally, along the z-axis. I recommend that the comments be read from inside
    to outside. */
    return std::lerp(
        /* Perform BIlinear interpolation between the four corners on the bottom face
        of the grid cube (the four corners with the minimum z-coordinate). This is
        equivalent to first interpolating along the x-axis and then interpolating along
        the y-axis (or first along the y-axis and then along the x-axis; both will give
        an equivalent result). I choose to interpolate along the x-axis first.
        
        Overall, this bilinear interpolation effectively computes the Value Noise at the point
        `min_coordinate_corner + GRID_CELL_SIZE * (interpolation_parameters.x,
        interpolation_parameters.y, 0)`. */
        std::lerp(
            /* Perform linear interpolation between the first adjacent pair of corners along the
            x-axis. Because the interpolation is along the x-axis, the interpolation parameter is
            `interpolation_parameters.x`. This effectively computes the Value Noise at the point
             `min_coordinate_corner + GRID_CELL_SIZE * (interpolation_parameters.x, 0, 0)`. */
            std::lerp(cell_corner_values[0], cell_corner_values[1], interpolation_parameters.x),
            /* Perform linear interpolation between the second adjacent pair of corners along the
            x-axis. Because the interpolation is along the x-axis, the interpolation parameter is
            `interpolation_parameters.x`. This effectively computes the Value Noise at the point
             `min_coordinate_corner + GRID_CELL_SIZE * (interpolation_parameters.x, 1, 0)`. */
            std::lerp(cell_corner_values[2], cell_corner_values[3], interpolation_parameters.x),
            /* Now, we take the above two interpolated values, and interpolate between them along
            the y-axis. Because the interpolation is performed along the y-axis, the interpolation
            parameter is `interpolation_parameters.y`. */
            interpolation_parameters.y
        ),
        /* Perform BIlinear interpolation between the four corners on the bottom face
        of the grid cube (the four corners with the minimum z-coordinate). This is
        equivalent to first interpolating along the x-axis and then interpolating along
        the y-axis (or first along the y-axis and then along the x-axis; both will give
        an equivalent result). Again, I choose to interpolate along the x-axis first.
        
        Overall, this bilinear interpolation effectively computes the Value Noise at the point
        `min_coordinate_corner + GRID_CELL_SIZE * (interpolation_parameters.x,
        interpolation_parameters.y, 1)`. */
        std::lerp(
            /* Perform linear interpolation between the third adjacent pair of corners along the
            x-axis. Because the interpolation is along the x-axis, the interpolation parameter is
            `interpolation_parameters.x`. This effectively computes the Value Noise at the point
             `min_coordinate_corner + GRID_CELL_SIZE * (interpolation_parameters.x, 0, 1)`. */
            std::lerp(cell_corner_values[4], cell_corner_values[5], interpolation_parameters.x),
            /* Perform linear interpolation between the first adjacent pair of corners along the
            x-axis. Because the interpolation is along the x-axis, the interpolation parameter is
            `interpolation_parameters.x`. This effectively computes the Value Noise at the point
             `min_coordinate_corner + GRID_CELL_SIZE * (interpolation_parameters.x, 1, 1)`. */
            std::lerp(cell_corner_values[6], cell_corner_values[7], interpolation_parameters.x),
            /* Now, we take the above two interpolated values, and interpolate between them along
            the y-axis. Because the interpolation is performed along the y-axis, the interpolation
            parameter is `interpolation_parameters.y`. */
            interpolation_parameters.y
        ),
        /* Finally, take the two bilinearly interpolated values from above, and perform a final
        linear interpolation of them along the z-axis. The interpolation parameter is given by
        `interpolation_parameters.z`, and the result of this linear interpolation is the Value
        Noise at the given point `p`, so we return it. */
        interpolation_parameters.z
    );
}

/* Returns the result of applying a rigid transformation to the point `p`. */
auto get_transformed(const Point3D &p) {
    /* The rigid transformation is represented by a 3x3 matrix, because we are in 3D space.
    This is the same 3x3 matrix that ssloy uses in their tutorial. The matrix which `p` is
    multiplied by is an orthogonal matrix; its rows and columns form an orthonormal basis.
    
    I believe this means it is a rigid transformation (aka an isometry), which is guaranteed to be
    a rotation, reflection, or a rotation + a reflection about a plane normal to the that is being
    applied to `p`.
    
    My intuitive understanding of this is that it essentially reinterprets the components of `p` as
    if they were the coordinates in a different basis of 3D space. That is, if p = (x, y, z) and the
    column vectors were v1, v2, v3, then we would be returning x * v1 + y * v2 + z * v3. */
    return Vec3D{
        dot(p, Vec3D{0.00,  0.80,  0.60}),
        dot(p, Vec3D{-0.80,  0.36, -0.48}),
        dot(p, Vec3D{-0.60, -0.48,  0.64})
    };
}

/* Given a point `p`, returns the corresponding Fractal Brownian Motion noise value. */
auto fractal_brownian_motion(const Point3D &p) {
    /* Here, we will not actually return the fractal noise at the point `p` itself. Instead, we will
    first apply a rigid transformation to `p`, and then return the fractal noise at that transformed
    point. ssloy does this, but doesn't explain why. GPT-4 says this could help to "break grid
    alignment", a common problem where noise functions generate structures that look aligned to the
    coordinate axes.
    
    In actuality, I can't visually tell the difference between explosions generated with and without
    transforming `p` before passing it to `value_noise_at`. I'm also doubtful of GPT-4's answer,
    because the transformation that we apply consists of rotations and/or reflections, meaning that
    if the original structure was grid-aligned, then the resulting structure should also be aligned,
    just along the direction of the basis vectors rather than the x/y/z-axes, right? However, I'm
    keeping the transformation in, for the sole reason that removing the transformation changes the
    explosion (because it means we are passing in different inputs to `value_noise_at`), and I don't
    want to spend more time looking for a new random seed I like. */
    auto transformed_point = get_transformed(p);

    auto ret = 0.;
    /* We add rescaled Value Noise to itself to produce what is called Fractal Noise.
    At each step, we halve the amplitude (the number we multiply by on the left), and
    we scale the frequency (which means we multiply `transformed_point`) by some
    number.
    
    Note that it is common practice to always double the frequency, but we do not
    do that here. This is because always doubling the frequency seems to result
    in explosions with very smooth surfaces, so much so that they almost look like
    blobs (likely because the consistent scaling of frequency causes peaks and
    troughs to line up somehow). Meanwhile, using different scaling factors disrupts
    this, resulting in more complex and irregular surfaces that more accurately resemble
    an explosion cloud. */
    ret += 0.5 * value_noise_at(transformed_point);
    ret += 0.25 * value_noise_at(transformed_point *= 2.32);
    ret += 0.125 * value_noise_at(transformed_point *= 3.03);
    ret += 0.0625 * value_noise_at(transformed_point *= 2.61);
    /* ^^ Four iterations gives a good enough level of detail. If you want more detail,
    add more iterations (and you may also want to increase the image resolution, or
    else the higher level of detail will not be discernable). */

    /* After the above, `ret` is in the range [0, 1 - 0.0625] (because the return value of
    `value_noise_at` is in the range [0, 1]), so we scale it to the range [0, 1]. This
    is equivalent to dividing `ret` by (1 - 0.0625). */
    return ret / (1 - 0.0625);
}

#endif