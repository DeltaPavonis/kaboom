#include <fstream>
#include <cstddef>
#include <vector>
#include <mutex>
#include "fractal_brownian_motion.h"

/* Dimensions of the rendered image */
constexpr size_t IMAGE_ROWS = 1080, IMAGE_COLS = 1080;

/* The center and radius of the implicit sphere. */
constexpr Point3D SPHERE_CENTER{0, 0, 0};
double SPHERE_RADIUS = 1.5;
/* Our noise values will fall in the range [0, `NOISE_AMPLITUDE`].
(The noise function always generates values in [0, 1], which we scale to [0, `NOISE_AMPLITUDE`]
by multiplying the returned noise values by `NOISE_AMPLITUDE`; see
`calculate_displaced_radius_towards`). */
constexpr double NOISE_AMPLITUDE = 1.;

/* Position of the single point light */
constexpr Point3D POINT_LIGHT_POS{10, 10, 10};

/* Properties of the camera */
constexpr Point3D CAMERA_CENTER{0, 0, 3};
constexpr double VERTICAL_FOV_RADIANS = std::numbers::pi / 3;  /* 60 degrees */
constexpr Vec3D BACKGROUND_COLOR{0.2, 0.7, 0.8};

/* This function introduces a displacement map on our sphere. Given any point `p`,
it returns the DISPLACED radius of the sphere at that point. Just as a point
`p` lies inside a regular sphere iff the distance from the sphere's center to `p`
is at most the radius of that sphere, a point `p` lies inside the DISPLACED sphere
iff the distance from the sphere's center to `p` is at most the DISPLACED radius
of the sphere towards that point `p`. */
auto calculate_displaced_radius_towards(const Point3D &p) {

    /* Previously, we first projected `p` to a point on the surface of the sphere, and only
    then did we apply our displacement map to that point on the surface of the sphere. By
    restricting the domain of the displacement map to points on the surface of the sphere,
    we ensured that the displacement map could only affect the surface of the sphere. And
    because the displacement map was continuous, that meant the displaced surface of the
    sphere was guaranteed to be continuous as well. In general, first projecting `p` to the
    surface of the sphere is suitable for creating bumps/indentations and other modifications
    to the sphere's surface that are still attached to the original sphere itself.
    
    However, this is not desirable for an explosion, which may have multiple small "clouds"
    that are detached from the original explosion cloud. As a result, we will now no longer
    project `p` onto the surface of the sphere before applying the displacement map to it.
    Instead, we will apply the displacement map to `p` directly. This allows the map to
    generate disconnected components and structures from the original sphere, which
    is what we want. */

    /* Now, our displacement map is defined in terms of Fractal Brownian Motion; that is,
    in terms of Fractal Noise. Specifically, at a point `p`, we will displace INWARD from
    the surface by `fractal_brownian_motion(p * 3.4) * NOISE_AMPLITUDE`. By the properties
    of Fractal Brownian Motion, this will result in a natural-looking cloud shape, which
    is exactly what we want.
    
    Note that the displacement depends on `fractal_brownian_motion(p * 3.4)`, rather than
    just `fractal_brownian_motion(p)`. This is what ssloy does, and I believe it is because
    Fractal Brownian Motion (and noise in general) generates similar values at points very
    close together. So, by multiplying `p` by 3.4, we increase the difference between the
    inputs to `fractal_brownian_motion` for points close together, resulting in less "flat"
    and more complex surfaces. Try adjusting the constant 3.4 to see for yourself; if it
    is too low, the explosion becomes just one big blob, and if it is too large, the
    explosion cloud becomes a disorganized mess (because we still WANT to have the
    property where points close together have similar displacements. If you multiply
    `p` by a constant that is too large, then even for points that are close together,
    we will be giving widely varying inputs to `fractal_brownian_motion`, resulting in a
    cloud with no inherent "connection"; the explosion will look like a collection of very
    tiny clouds, each one representing a group of points very near to each other; so near
    that even with the large multiplier, they were generated in the same connected component
    by the Fractal Brownian Motion). */
    auto val = fractal_brownian_motion(p * 3.4);  /* See above comment about the 3.4 */
    auto displacement = -val * NOISE_AMPLITUDE;  /* Displacement is in [-NOISE_AMPLITUDE, 0] */

    /* Now, the displacement map maps arbitrary points in space to their displacements.
    As before, we return `SPHERE_RADIUS` + `displacement` for the displaced radius of the
    displacement-mapped sphere towards the point `p`. */
    return SPHERE_RADIUS + displacement;
}

/* Returns the signed distance from the point `p` to the surface of the sphere.
By "signed distance", we mean the value whose magnitude equals the shortest distance
from `p` to the sphere's surface, and whose sign specifies whether or not `p` is
inside the sphere. Specifically, if the signed distance is negative, positive, or
zero, then `p` lies inside, outside, and on the surface of the sphere, respectively.
Thus, the point `p` lies on or inside the sphere if and only if the signed distance
from `p` to the sphere is nonpositive. */
auto signed_distance_from_sphere(const Point3D &p) {
    /* The shortest distance from a point `p` to a sphere with center `C` and radius `r`
    is simply |dist(C, p) - r| (because it just equals the difference between the distance
    from the center to `p`, and the distance from the center to the sphere's surface (which
    is just the radius `r`)). As a result, the SIGNED distance is just dist(C, p) - r, because
    its magnitude is |dist(C, p) - r|, which we just established was the shortest distance
    from `p` to the sphere's surface, and because the sign of dist(C, p) - r tells us when
    `p` is inside, outside, or on the sphere (dist(C, p) - r is negative iff dist(C, p) < r
    iff p is inside the sphere, dist(C, p) - r is positive iff dist(C, p) > r iff p is outside
    the sphere, and dist(C, p) - r is zero iff dist(C, p) = r iff p is on the sphere's surface,
    as desired).

    Finally, because our sphere is centered at the origin, dist(C, p) is just equal to the
    magnitude of p. Thus, the signed distance from `p` to the sphere is given by
    p.mag() - SPHERE_RADIUS, which is what we will return.
    
    After adding displacement mapping, the only modification needed to return the signed distance
    to the displacement-mapped sphere is simply to use the displaced radius towards `p` in place
    of the original radius `SPHERE_RADIUS`.*/
    return (p - SPHERE_CENTER).mag() - calculate_displaced_radius_towards(p);
}

/* Returns the closest hit point of the ray with origin `ray_origin` and direction `ray_dir`
with the sphere we are rendering, which is centered at the origin, and has radius given by
`SPHERE_RADIUS`. */
std::optional<Point3D> closest_sphere_hit_point(const Point3D &ray_origin, const Vec3D &ray_dir) {

    /* Because our displacement map never increases the radius of the sphere towards any point
    (as the `displacement` variable in `calculate_displaced_radius_towards` is always negative),
    we know that a ray intersects with our explosion only if it intersects with the original
    non-displacement-mapped sphere, which is the sphere centered at SPHERE_CENTER and with
    radius SPHERE_RADIUS. Thus, for optimization purposes, we will first do some basic 3D
    geometry and math to check if the ray intersects with the original sphere; if it does not,
    then we immediately know that it will not intersect the explosion, and so we will return
    an empty `std::optional<Point3D>` as a result. */
    {
        /* Set up quadratic formula calculation */
        auto center_to_origin = ray_origin - SPHERE_CENTER;
        auto a = dot(ray_dir, ray_dir);
        auto b_half = dot(ray_dir, center_to_origin);
        auto c = dot(center_to_origin, center_to_origin) - SPHERE_RADIUS * SPHERE_RADIUS;
        auto discriminant_quarter = b_half * b_half - a * c;

        /* Quadratic has no solutions whenever the discriminant is negative; this means
        the given ray does not intersect the sphere */
        if (discriminant_quarter < 0) {return {};}
    }

    /* `ray_unit_dir` = An UNIT VECTOR parallel to the given ray's direction. */
    auto ray_unit_dir = ray_dir.unit_vector();

    /* To compute the closest hit point, we could find the hit point analytically by considering
    the equation of a sphere in 3D space and then solving for the hit times. A less accurate but
    simpler approach is to simply test increasing distances from the origin along our given ray
    until we first hit or go inside the sphere; when that happens, we just return the current
    point, and that will be the closest hit point to the sphere.
    
    Here, we test distances starting from 0 (obviously) in increments of 0.01. In our scene,
    every closest hit point has distance less than 5 from the camera center, so we only need
    to test distances less than 5. */
    for (auto distance_along_ray = 0.; distance_along_ray < 5; distance_along_ray += 0.05) {

        /* Check if the point at a distance of `distance_along_ray` along the given ray is
        on or inside the sphere. If it does hit the sphere, then it is the closest hit point
        of this ray with the sphere, so we return it immediately. */
        auto curr_point = ray_origin + distance_along_ray * ray_unit_dir;
        
        /* `curr_point` is on or inside the sphere iff its signed distance to the sphere
        is nonpositive, by the definition of signed distance (see the comments for
        `signed_distance_from_sphere()`). */
        if (signed_distance_from_sphere(curr_point) <= 0) {
            return curr_point;
        }
    }

    /* If the ray never hits the sphere, we return an empty `std::optional` object. */
    return {};
}

/* Approximates the unit surface normal at `hit_point` on the displacement-mapped sphere
by using the method of finite differences. */
auto get_unit_surface_normal(const Point3D &hit_point) {
    /* The key idea here is that finding the (outward) surface normal at a point `hit_point`
    is equivalent to finding the gradient of the signed distance function at that same point.
    This is because the gradient points in the direction of steepest increase for a function,
    which, for the signed distance function (and for distance functions in general), is
    outward and perpendicular to the surface. Thus, it suffices to find the gradient of
    the signed distance function at `hit_point`, and then normalize it. And to do this, we
    need to find the partial derivatives of the signed distance function with respect to the
    x-, y-, and z- coordinates, at the point `hit_point`.
 
    We approximate the partial derivatives of the signed distance with respect to the x-, y-
    and z-coordinates at the point `hit_point` by using the method of finite differences. The
    most basic form of the method of finite differences, which we use here, is very simple:
    it just states that to approximate f'(x), you evaluate (f(x + h) - f(x)) / h for some
    very small positive real number h. This works because the definition of the derivative
    f'(x) is just the limit, as h approaches 0, of (f(x + h) - f(x)) / h.
    
    Using this idea, we approximate the partial derivative of the signed distance with respect
    to some axis (say the x-axis) at the point `hit_point` by evaluating the expression
    signed_distance(hit_point shifted by some epsilon along the x-axis) - signed_distance(hit_point)
    While we technically need to divide this quantity by `epsilon` in order for it to approximate
    the partial derivative, we don't need that here. This is because again, we just need to return
    the normalized gradient; all that matters is that the ratio between the partial derivatives
    is preserved, and so we can just leave out dividing by epsilon for all three of the partial
    derivatives, because not dividing does not change the ratio between the partial derivatives. */

    constexpr double epsilon = 0.1;
    auto initial_dist = signed_distance_from_sphere(hit_point);
    auto partial_x = signed_distance_from_sphere(hit_point + Vec3D{epsilon, 0, 0}) - initial_dist;
    auto partial_y = signed_distance_from_sphere(hit_point + Vec3D{0, epsilon, 0}) - initial_dist;
    auto partial_z = signed_distance_from_sphere(hit_point + Vec3D{0, 0, epsilon}) - initial_dist;

    /* Return the normalized gradient at `hit_point`, which is just the unit vector of
    {partial_x, partial_y, partial_z}. This is our unit approximated surface normal at the point
    `hit_point`. */
    return Vec3D{partial_x, partial_y, partial_z}.unit_vector();
}

int main()
{
    /* The seed used for the images that come with this repository. Delete or modify this
    line if you want to set your own seed, or if you just want the program to use a
    random seed (the program automatically generates a random seed if no seed is explicitly
    provided by the user; see the next line of code). */
    seed = 0x1.167d98ef1184p+19;

    /* Generate a random seed if a seed was not provided by the user */
    generate_seed_if_not_set();

    /* Print the seed in hexadecimal floating-point format to ensure exactness */
    std::cout << std::showbase << std::hexfloat << "Seed: " << *seed << std::endl;
    std::cout << std::noshowbase << std::defaultfloat;  /* Undo `std::showbase`/`std::hexfloat` */
 
    /* Render image of dimensions IMAGE_COLS x IMAGE_ROWS */
    std::vector image(IMAGE_ROWS, std::vector<Vec3D>(IMAGE_COLS));

    /* Now, let's set up the viewport. For maximum simplicity, we will just set the viewport
    height and width to be equal to `IMAGE_ROWS` and `IMAGE_COLS`, respectively. */
    const double VIEWPORT_HEIGHT = IMAGE_ROWS, VIEWPORT_WIDTH = IMAGE_COLS;
    /* Every other viewport property - the focal length (shortest distance from the camera
    center to the viewport, the top-left corner of the viewport, the delta for moving one
    row or one column, and the location of the top-left pixel - are all determined by
    the viewport height, viewport width, camera center, and camera field-of-view (vertical FOV
    here; you can also use horizontal FOV). */
    /* `FOCAL_LENGTH` is calculated by the usual formula, as is `VIEWPORT_TOP_LEFT_CORNER`. */
    const auto FOCAL_LENGTH = VIEWPORT_HEIGHT / (2. * std::tan(VERTICAL_FOV_RADIANS / 2.));
    const auto VIEWPORT_TOP_LEFT_CORNER = CAMERA_CENTER +
        Point3D{
            -VIEWPORT_WIDTH / 2.,  /* Need to move left from the camera center, so negative x */
            VIEWPORT_HEIGHT / 2,   /* Need to move up from the camera center, so positive y */
            -FOCAL_LENGTH  /* By right-handed coordinates, camera points towards negative z */
        };
    /* Because we set the viewport height/width equal to the image height/width, the row and
    column deltas are just {0, -1, 0} (negative 1 because we are moving down a row, and so the
    y-coordinate decreases) and {1, 0, 0}, respectively. */
    const Point3D PIXEL_ROW_DELTA{0, -1, 0}, PIXEL_COL_DELTA{1, 0, 0};
    /* The center of the top-left pixel's square is calculated by the usual formula */
    const auto PIXEL00_LOC = VIEWPORT_TOP_LEFT_CORNER + PIXEL_ROW_DELTA / 2 + PIXEL_COL_DELTA / 2;

    /* Variables needed to continuously update and print out the number of rows left in the
    render loop below */
    auto rows_left = IMAGE_ROWS;
    std::mutex update_rows_left_mtx;
    
    /* Render Loop */
    #pragma omp parallel for schedule(dynamic, 8)  /* Parallelize across rows */
    for (size_t row = 0; row < IMAGE_ROWS; ++row) {
        for (size_t col = 0; col < IMAGE_COLS; ++col) {

            /* To compute the color for the current pixel, we need to shoot a ray originating from
            CAMERA_CENTER, and through the center of the square on the viewport corresponding to
            the current pixel. 
            
            Here, `ray_endpoint` equals the point at the center of the current pixel's square on
            the viewport. The ray direction is then found by taking the difference between the
            ray's endpoint and its origin (which is `CAMERA_CENTER`).
            
            You can observe that all this is just equivalent to setting `ray_dir` to
            Vec3D {
                i - VIEWPORT_WIDTH / 2 + 0.5 = i - IMAGE_WIDTH / 2 + 0.5,
                VIEWPORT_HEIGHT / 2 - 0.5 - j = IMAGE_HEIGHT / 2 - 0.5 - j,
                FOCAL_LENGTH = -VIEWPORT_HEIGHT / 2 * tan(VERTICAL_FOV_RADIANS / 2)
                             = -(double) IMAGE_HEIGHT / 2 * tan(VERTICAL_FOV_RADIANS / 2)
            }
            which is what they do in the tutorial. By understanding how they set up the viewport,
            which we did before this `for`-loop, you can understand why that code works. */
            auto ray_endpoint = PIXEL00_LOC
                                + PIXEL_ROW_DELTA * static_cast<double>(row)
                                + PIXEL_COL_DELTA * static_cast<double>(col);
            auto ray_dir = ray_endpoint - CAMERA_CENTER;

            /* Shoot a camera ray with origin `CAMERA_CENTER` and direction given by `ray_dir`;
            if it hits the sphere, then the current pixel will be white, and otherwise it will
            be the default background color (a blueish shade). */
            if (auto hit_point = closest_sphere_hit_point(CAMERA_CENTER, ray_dir); hit_point) {
                /* This camera ray hit the sphere */

                /* Implement Lambertian (diffuse) reflectance for the sphere; the brightness
                at a certain point p on the sphere is proportional to the cosine of the angle
                between the surface normal at that hit point p, and the vector from p towards
                the point light. Note that this is different from Lambert's Cosine Law (I think),
                which uses the cosine of the angle between the surface normal and the vector
                from the observer towards the hit point. According to GPT-4, the latter is used
                for physically modeling surface emission and scattering of light, while the former
                models reflection of light; these are different interpretations of Lambert's
                Cosine Law, and the former is used for computer graphics. TODO: Fact-check this.
                
                To compute the cosine of the angle between the surface normal and the direction
                towards the point light, we observe that for two vectors a, b, we have
                dot(a, b) = |a||b|cos(theta). Thus, if |a| = |b| = 1, then dot(a, b) = cos(theta),
                which is what we want to compute. As a result, here, cos(theta) is equal to the
                dot product of the UNIT surface normal at the hit point, and the UNIT direction
                from the hit point towards the light. */

                /* Because our object is no longer a sphere (due to our displacement mapping),
                the unit surface normal at `hit_point` can no longer be easily directly
                computed. As a result, we introduce `get_unit_surface_normal`, which numerically
                approximates the unit surface normal at `hit_point` (see the comments for it
                above). For our purposes, this is good enough. */
                auto unit_surface_normal = get_unit_surface_normal(*hit_point);
                auto unit_dir_to_light = (POINT_LIGHT_POS - *hit_point).unit_vector();
                auto brightness_factor = std::clamp(
                    /* Calculate the cosine of the angle between the surface normal at the hit
                    point, and the direction from the hit point to the light */
                    dot(unit_surface_normal, unit_dir_to_light),
                    /* Then clamp that value, because it may be negative or exceed 1. The
                    tutorial chooses to set a lower bound of 0.4, so I will do that as well. */
                    0.4, 1.
                );

                /* After applying the displacement mapping, we don't need our spatial texture
                anymore. We just let Lambertian reflectance from the point light handle the
                colors. */

                /* We say that the intrinsic color of every point on our sphere is just pure white
                (Vec3D{1, 1, 1}). This color is scaled by the `brightness_factor`, which is
                determined by our point light and by Lambertian reflectance above. */
                image[row][col] = Vec3D{1, 1, 1} * brightness_factor;
            } else {
                /* This camera ray did not hit the sphere, so this pixel's color will be equal
                to the background color */
                image[row][col] = BACKGROUND_COLOR;
            }
        }

        /* Update and then print the number of rows left to render */
        std::lock_guard guard(update_rows_left_mtx);
        std::cout << "\r" << --rows_left << " rows left" << ' ' << std::flush;
        /* ^^ Note that there is a space at the end of the line; this is intentional and necessary.
        When the length of `rows_left` decreases by 1, then the space at the end of the line serves
        to clear the last character of the line, so that it doesn't show up later (that is, if
        we did not print a space at the end of the line, then if we started with "1000 rows left"
        the next iteration would print "999 rows leftt" because the previous last character "t"
        would still remain on the line). */
    }

    /* Output the image as a PPM file to the location specified by `output_file_path` */
    std::string output_file_path = "rendered_image.ppm";
    std::ofstream fout(output_file_path);
    fout << "P3\n" << IMAGE_COLS << " " << IMAGE_ROWS << "\n255\n";
    for (size_t row = 0; row < IMAGE_ROWS; ++row) {
        for (size_t col = 0; col < IMAGE_COLS; ++col) {
            /* Scale each real component in [0, 1] to an integer component in [0, 255] */
            int r = static_cast<int>(image[row][col].x * 255.99999);
            int g = static_cast<int>(image[row][col].y * 255.99999);
            int b = static_cast<int>(image[row][col].z * 255.99999);
            fout << r << " " << g << " " << b << '\n';
        }
    }

    std::cout << "\rSuccessfully outputted image to " << output_file_path << std::endl;
    
    return 0;
}