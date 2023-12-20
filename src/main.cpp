#include <fstream>
#include <cstddef>
#include <vector>
#include <numbers>
#include <optional>
#include "vec3d.h"

/* Dimensions of the rendered image */
constexpr size_t IMAGE_ROWS = 1000, IMAGE_COLS = 1000;

/* Radius of the single sphere, which is centered at the origin */
constexpr double SPHERE_RADIUS = 1.5;

/* Properties of the camera */
constexpr Point3D CAMERA_CENTER{0, 0, 3};
constexpr double VERTICAL_FOV_RADIANS = std::numbers::pi / 3;
constexpr Vec3D BACKGROUND_COLOR{0.2, 0.7, 0.8};

/* Returns the closest hit point of the ray with origin `ray_origin` and direction `ray_dir`
with the sphere we are rendering, which is centered at the origin, and has radius given by
`SPHERE_RADIUS`. */
std::optional<Point3D> closest_sphere_hit_point(const Point3D &ray_origin, const Vec3D &ray_dir) {

    /* `ray_unit_dir` = An UNIT VECTOR parallel to the given ray's direction. */
    auto ray_unit_dir = ray_dir.unit_vector();

    /* To compute the closest hit point, we could find the hit point analytically by considering
    the equation of a sphere in 3D space and then solving for the hit times. A less accurate but
    simpler approach is to simply test increasing distances from the origin along our given ray
    until we first hit or go inside the sphere; when that happens, we just return the current
    point, and that will be the closest hit point to the sphere.
    
    Here, we test distances starting from 0 (obviously) in increments of 0.01. */
    for (double distance_along_ray = 0; distance_along_ray < 100; distance_along_ray += 0.01) {

        /* Check if the point at a distance of `distance_along_ray` along the given ray is
        on or inside the sphere. To do this, observe that a point (x, y, z) is on/inside the
        sphere if and only if x^2 + y^2 + z^2 <= SPHERE_RADIUS^2, because the sphere is
        centered at the origin. If the current point (which is given by taking ray_origin,
        then adding distance_along_ray times the UNIT ray direction) does hit the sphere,
        then it is the closest hit point of this ray with the sphere, so we just
        return it immediately. */
        auto curr_point = ray_origin + distance_along_ray * ray_unit_dir;
        if (curr_point.mag_squared() <= SPHERE_RADIUS * SPHERE_RADIUS) {
            return curr_point;
        }
    }

    /* If the ray never hits the sphere, we return an empty `std::optional` object. */
    return {};
}

int main()
{
    /* Render image of dimensions IMAGE_ROWS x IMAGE_COLS */
    std::vector image(IMAGE_ROWS, std::vector<Vec3D>(IMAGE_COLS));

    /* Now, let's set up the viewport. For maximum simplicity, we will just set the viewport
    height and width to be equal to `IMAGE_ROWS` and `IMAGE_COLS`, respectively. */
    const double VIEWPORT_HEIGHT = IMAGE_ROWS, VIEWPORT_WIDTH = IMAGE_COLS;
    /* Every other viewport property - the focal length (shortest distance from the camera
    center to the viewport, the top-left corner of the viewport, the delta for moving one
    row or one column, and the location of the top-left pixel - are all determined by
    the viewport height, viewport width, and camera center. */
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

    /* Render loop */
    #pragma omp parallel for schedule(dynamic)
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
                image[row][col] = Vec3D{1, 1, 1};  /* The sphere's color is pure white */
            } else {
                /* This camera ray did not hit the sphere, so this pixel's color will be equal
                to the background color */
                image[row][col] = BACKGROUND_COLOR;
            }
        }
    }

    /* Output the image as a PPM file called "image.ppm" */
    std::ofstream fout("rendered_image.ppm");
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

    return 0;
}