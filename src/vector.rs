// Copyright (c) 2024 Ken Barker

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

//! The `vector` module contains functions for performing great circle
//! calculations using vectors to represent points and great circle poles on a
//! unit sphere.

use crate::great_circle;
use crate::Vector3d;
use angle_sc::trig;
use angle_sc::{Angle, Radians};

/// The minimum length of a vector to normalize.
pub const MIN_LENGTH: f64 = 16384.0 * core::f64::EPSILON;

/// The minimum norm of a vector to normalize.
pub const MIN_NORM: f64 = MIN_LENGTH * MIN_LENGTH;

/// Determine whether a `Vector3d` is a unit vector.
///
/// returns true if `Vector3d` is a unit vector, false otherwise.
#[must_use]
pub fn is_unit(a: &Vector3d) -> bool {
    const MIN_POINT_SQ_LENGTH: f64 = 1.0 - 12.0 * core::f64::EPSILON;
    const MAX_POINT_SQ_LENGTH: f64 = 1.0 + 12.0 * core::f64::EPSILON;

    (MIN_POINT_SQ_LENGTH..=MAX_POINT_SQ_LENGTH).contains(&(a.norm()))
}

/// Calculate the square of the Euclidean distance between two Points.
/// Note: points do NOT need to be valid Points.  
/// @post for unit vectors: result <= 4
#[must_use]
pub fn sq_distance(a: &Vector3d, b: &Vector3d) -> f64 {
    (b - a).norm_squared()
}

/// Calculate the shortest (Euclidean) distance between two Points.  
/// @post for unit vectors: result <= 2
#[must_use]
pub fn distance(a: &Vector3d, b: &Vector3d) -> f64 {
    (b - a).norm()
}

/// Determine whether two `Vector3d`s are orthogonal (perpendicular).
///
/// returns true if a and b are orthogonal, false otherwise.
#[must_use]
pub fn are_orthogonal(a: &Vector3d, b: &Vector3d) -> bool {
    const MAX_LENGTH: f64 = 4.0 * core::f64::EPSILON;

    (-MAX_LENGTH..=MAX_LENGTH).contains(&(a.dot(b)))
}

/// Determine whether point a is West of point b.  
/// It calculates and compares the perp product of the two points.
/// * `a`, `b` - the points.
///
/// returns true if a is West of b, false otherwise.
#[must_use]
pub fn is_west_of(a: &Vector3d, b: &Vector3d) -> bool {
    // Compare with -epsilon to handle floating point errors
    b.xy().perp(&a.xy()) <= -core::f64::EPSILON
}

/// Calculate the right hand pole vector of a Great Circle from an initial
/// position and an azimuth.  
/// See: <http://www.movable-type.co.uk/scripts/latlong-vectors.html#Great%20circles>
/// * `lat` - start point Latitude.
/// * `lon` - start point Longitude.
/// * `azi` - start point azimuth.
///
/// returns the right hand pole vector of the great circle.
#[must_use]
pub fn calculate_pole(lat: Angle, lon: Angle, azi: Angle) -> Vector3d {
    let x = trig::UnitNegRange::clamp(
        lon.sin().0 * azi.cos().0 - lat.sin().0 * lon.cos().0 * azi.sin().0,
    );
    let y = trig::UnitNegRange::clamp(
        0.0 - lon.cos().0 * azi.cos().0 - lat.sin().0 * lon.sin().0 * azi.sin().0,
    );
    let z = trig::UnitNegRange(lat.cos().0 * azi.sin().0);

    Vector3d::new(x.0, y.0, z.0)
}

/// Calculate the azimuth at a point on the Great Circle defined by pole.
/// * `point` - the point.
/// * `pole` - the right hand pole of the Great Circle.
///
/// returns the azimuth at the point on the great circle.
#[must_use]
pub fn calculate_azimuth(point: &Vector3d, pole: &Vector3d) -> Angle {
    const MAX_LAT: f64 = 1.0 - 2.0 * core::f64::EPSILON;

    let sin_lat = point.z;
    // if the point is close to the North or South poles, azimuth is 180 or 0.
    if MAX_LAT <= libm::fabs(sin_lat) {
        return if 0.0 < sin_lat {
            Angle::from_y_x(0.0, -1.0)
        } else {
            Angle::default()
        };
    }

    Angle::from_y_x(pole.z, pole.xy().perp(&point.xy()))
}

/// Calculate the direction vector along a Great Circle from an initial
/// position and an azimuth.  
/// See: Panou and Korakitis equations: 30, 31, & 32a
/// <https://arxiv.org/abs/1811.03513>
/// * `lat` - start point Latitude.
/// * `lon` - start point Longitude.
/// * `azi` - start point azimuth.
///
/// returns the direction vector at the point on the great circle.
#[must_use]
pub fn calculate_direction(lat: Angle, lon: Angle, azi: Angle) -> Vector3d {
    let x = trig::UnitNegRange::clamp(
        0.0 - lat.sin().0 * lon.cos().0 * azi.cos().0 - lon.sin().0 * azi.sin().0,
    );
    let y = trig::UnitNegRange::clamp(
        0.0 - lat.sin().0 * lon.sin().0 * azi.cos().0 + lon.cos().0 * azi.sin().0,
    );
    let z = trig::UnitNegRange(lat.cos().0 * azi.cos().0);

    Vector3d::new(x.0, y.0, z.0)
}

/// Calculate the direction vector of a Great Circle arc.
/// * `a` - the start point.
/// * `pole` - the pole of a Great Circle.
///
/// returns the direction vector at the point on the great circle.
#[must_use]
pub fn direction(a: &Vector3d, pole: &Vector3d) -> Vector3d {
    pole.cross(a)
}

/// Calculate the position of a point along a Great Circle arc.
/// * `a` - the start point.
/// * `dir` - the direction vector of a Great Circle at a.
/// * `distance` - the a Great Circle as an Angle.
///
/// returns the position vector at the point on the great circle.
#[must_use]
pub fn position(a: &Vector3d, dir: &Vector3d, distance: Angle) -> Vector3d {
    distance.cos().0 * a + distance.sin().0 * dir
}

/// Calculate the direction vector of a Great Circle rotated by angle.
/// * `dir` - the direction vector of a Great Circle arc.
/// * `pole` - the pole of a Great Circle.
/// * `angle` - the angle to rotate the direction vector by.
///
/// returns the direction vector at the point on the great circle
/// rotated by angle.
#[must_use]
pub fn rotate(dir: &Vector3d, pole: &Vector3d, angle: Angle) -> Vector3d {
    position(dir, pole, angle)
}

/// Calculate the position of a point rotated by angle at radius.
/// * `a` - the start point.
/// * `pole` - the pole of a Great Circle.
/// * `angle` - the angle to rotate the direction vector by.
/// * `radius` - the radius from the start point.
///
/// returns the position vector at angle and radius from the start point.
#[must_use]
pub fn rotate_position(a: &Vector3d, pole: &Vector3d, angle: Angle, radius: Angle) -> Vector3d {
    position(a, &rotate(&direction(a, pole), pole, angle), radius)
}

/// The sine of the across track distance of a point relative to a Great Circle pole.  
/// It is simply the dot product of the pole and the point: pole . point
/// * `pole` - the Great Circle pole.
/// * `point` - the point.
///
/// returns the sine of the across track distance of point relative to the pole.
#[must_use]
fn sin_xtd(pole: &Vector3d, point: &Vector3d) -> trig::UnitNegRange {
    trig::UnitNegRange::clamp(pole.dot(point))
}

/// The across track distance of a point relative to a Great Circle pole.
/// * `pole` - the Great Circle pole.
/// * `point` - the point.
///
/// returns the across track distance of point relative to pole.
#[must_use]
pub fn cross_track_distance(pole: &Vector3d, point: &Vector3d) -> Radians {
    Radians(libm::asin(sin_xtd(pole, point).0))
}

/// The square of the Euclidean cross track distance of a point relative to a
/// Great Circle pole.
/// * `pole` - the Great Circle pole.
/// * `point` - the point.
///
/// returns the square of the euclidean distance of point relative to pole.
#[must_use]
pub fn sq_cross_track_distance(pole: &Vector3d, point: &Vector3d) -> f64 {
    let sin_d = libm::fabs(sin_xtd(pole, point).0);
    if sin_d < 2.0 * core::f64::EPSILON {
        0.0
    } else {
        2.0 * (1.0 - trig::swap_sin_cos(trig::UnitNegRange(sin_d)).0)
    }
}

/// Calculate the closest point on a plane to the given point.
/// * `pole` - the Great Circle pole (aka normal) of the plane.
/// * `point` - the point.
///
/// returns the closest point on a plane to the given point.
#[must_use]
fn calculate_point_on_plane(pole: &Vector3d, point: &Vector3d) -> Vector3d {
    let t = sin_xtd(pole, point);
    point - pole * t.0
}

/// Calculate the square of the Euclidean along track distance of a point
/// from the start of an Arc.
/// It is calculated using the closest point on the plane to the point.
/// * `a` - the start point of the Great Circle arc.
/// * `pole` - the pole of the Great Circle arc.
/// * `point` - the point.
/// returns the square of the Euclidean along track distance
#[must_use]
pub fn sq_along_track_distance(a: &Vector3d, pole: &Vector3d, point: &Vector3d) -> f64 {
    let plane_point = calculate_point_on_plane(pole, point);
    if plane_point.norm() < MIN_LENGTH {
        0.0
    } else {
        sq_distance(a, &(plane_point.normalize()))
    }
}

/// The sine of the along track distance of a point along a Great Circle arc.  
/// It is the triple product of the pole, a and the point:
/// (pole X a) . point = pole . (a X point)
/// * `a` - the start point of the Great Circle arc.
/// * `pole` - the pole of the Great Circle arc.
/// * `point` - the point.
///
/// returns the sine of the along track distance of point relative to the start
/// of a great circle arc.
#[must_use]
pub fn sin_atd(a: &Vector3d, pole: &Vector3d, point: &Vector3d) -> trig::UnitNegRange {
    trig::UnitNegRange::clamp(pole.cross(a).dot(point))
}

/// The square of `core::f64::EPSILON`.
pub const SQ_EPSILON: f64 = core::f64::EPSILON * core::f64::EPSILON;

/// The Great Circle distance of a point along the arc relative to a,
/// (+ve) ahead of a, (-ve) behind a.
/// * `a` - the start point of the Great Circle arc.
/// * `pole` - the pole of the Great Circle arc.
/// * `point` - the point.
///
/// returns the along track distance of point relative to the start of a great circle arc.
#[must_use]
pub fn along_track_distance(a: &Vector3d, pole: &Vector3d, point: &Vector3d) -> Radians {
    let sq_atd = sq_along_track_distance(a, pole, point);
    if sq_atd <= SQ_EPSILON {
        Radians(0.0)
    } else {
        Radians(libm::copysign(
            great_circle::e2gc_distance(libm::sqrt(sq_atd)).0,
            sin_atd(a, pole, point).0,
        ))
    }
}

/// Calculate Great Circle along and across track distances.
/// * `a` - the start point of the Great Circle arc.
/// * `pole` - the pole of the Great Circle arc.
/// * `p` - the point.
///
/// returns the along and across track distances of point relative to the
/// start of a great circle arc.
#[allow(clippy::similar_names)]
#[must_use]
pub fn calculate_atd_and_xtd(a: &Vector3d, pole: &Vector3d, p: &Vector3d) -> (Radians, Radians) {
    let mut atd = Radians(0.0);
    let mut xtd = Radians(0.0);

    let mut sq_d = sq_distance(a, p);
    // Point is close to start
    if 2.0 * SQ_EPSILON < sq_d {
        let sin_xtd = sin_xtd(pole, p);
        let abs_sin_xtd = libm::fabs(sin_xtd.0);

        if abs_sin_xtd < 1.0 - 2.0 * core::f64::EPSILON {
            // if the across track distance is significant
            if core::f64::EPSILON < abs_sin_xtd {
                // calculate the signed great circle across track distance
                xtd = Radians(libm::asin(sin_xtd.0));

                // the closest point on the plane of the pole to the point
                let plane_point = p - pole * sin_xtd.0;
                sq_d = if MIN_LENGTH < plane_point.norm() {
                    sq_distance(a, &(plane_point.normalize()))
                } else {
                    0.0
                };
            }
            atd = Radians(libm::copysign(
                great_circle::e2gc_distance(libm::sqrt(sq_d)).0,
                sin_atd(a, pole, p).0,
            ));
        } else {
            xtd = Radians(if 0.0 < sin_xtd.0 {
                core::f64::consts::FRAC_PI_2
            } else {
                -core::f64::consts::FRAC_PI_2
            });
        }
    }

    (atd, xtd)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::LatLong;
    use angle_sc::{is_within_tolerance, Degrees};

    #[test]
    fn test_point_lat_longs() {
        // Test South pole
        let lat_lon_south = LatLong::new(Degrees(-90.0), Degrees(180.0));
        let point_south = Vector3d::from(&lat_lon_south);

        let result = LatLong::from(&point_south);
        assert_eq!(-90.0, result.lat().0);
        // Note: longitude is now zero, since the poles do not have a Longitude
        assert_eq!(0.0, result.lon().0);

        // Test Greenwich equator
        let lat_lon_0_0 = LatLong::new(Degrees(0.0), Degrees(0.0));
        let point_0 = Vector3d::from(&lat_lon_0_0);
        assert!(is_unit(&point_0));
        assert_eq!(lat_lon_0_0, LatLong::from(&point_0));

        // Test IDL equator
        let lat_lon_0_180 = LatLong::new(Degrees(0.0), Degrees(180.0));
        let point_1 = Vector3d::from(&lat_lon_0_180);
        assert!(is_unit(&point_1));
        assert_eq!(lat_lon_0_180, LatLong::from(&point_1));
        assert_eq!(false, is_west_of(&point_0, &point_1));

        let lat_lon_0_m180 = LatLong::new(Degrees(0.0), Degrees(-180.0));
        let point_2 = Vector3d::from(&lat_lon_0_m180);
        assert!(is_unit(&point_2));
        // Converts back to +ve longitude
        assert_eq!(lat_lon_0_180, LatLong::from(&point_2));
        assert_eq!(false, is_west_of(&point_0, &point_2));

        let lat_lon_0_r3 = LatLong::new(Degrees(0.0), Degrees(3.0_f64.to_degrees()));
        let point_3 = Vector3d::from(&lat_lon_0_r3);
        assert!(is_unit(&point_3));
        let result = LatLong::from(&point_3);
        assert_eq!(0.0, result.lat().0);
        assert_eq!(3.0_f64.to_degrees(), result.lon().0);
        assert!(is_west_of(&point_0, &point_3));
        assert_eq!(false, is_west_of(&point_1, &point_3));

        let lat_lon_0_mr3 = LatLong::new(Degrees(0.0), Degrees(-3.0_f64.to_degrees()));
        let point_4 = Vector3d::from(&lat_lon_0_mr3);
        assert!(is_unit(&point_4));
        let result = LatLong::from(&point_4);
        assert_eq!(0.0, result.lat().0);
        assert_eq!(-3.0_f64.to_degrees(), result.lon().0);
        assert!(is_west_of(&point_1, &point_4));
    }

    #[test]
    fn test_point_distance() {
        let lat_lon_south = LatLong::new(Degrees(-90.0), Degrees(0.0));
        let south_pole = Vector3d::from(&lat_lon_south);

        let lat_lon_north = LatLong::new(Degrees(90.0), Degrees(0.0));
        let north_pole = Vector3d::from(&lat_lon_north);

        assert_eq!(0.0, sq_distance(&south_pole, &south_pole));
        assert_eq!(0.0, sq_distance(&north_pole, &north_pole));
        assert_eq!(4.0, sq_distance(&south_pole, &north_pole));

        assert_eq!(0.0, distance(&south_pole, &south_pole));
        assert_eq!(0.0, distance(&north_pole, &north_pole));
        assert_eq!(2.0, distance(&south_pole, &north_pole));

        // Greenwich equator
        let g_eq = Vector3d::new(1.0, 0.0, 0.0);

        // Test IDL equator
        let idl_eq = Vector3d::new(-1.0, 0.0, 0.0);

        assert_eq!(0.0, sq_distance(&g_eq, &g_eq));
        assert_eq!(0.0, sq_distance(&idl_eq, &idl_eq));
        assert_eq!(4.0, sq_distance(&g_eq, &idl_eq));

        assert_eq!(0.0, distance(&g_eq, &g_eq));
        assert_eq!(0.0, distance(&idl_eq, &idl_eq));
        assert_eq!(2.0, distance(&g_eq, &idl_eq));
    }

    #[test]
    fn test_calculate_azimuth_at_poles() {
        // Greenwich equator
        let g_eq = Vector3d::new(1.0, 0.0, 0.0);
        let south_pole = Vector3d::new(0.0, 0.0, -1.0);
        let result = calculate_azimuth(&south_pole, &g_eq);
        assert_eq!(Angle::default(), result);

        let north_pole = Vector3d::new(0.0, 0.0, 1.0);
        let result = calculate_azimuth(&north_pole, &g_eq);
        assert_eq!(Angle::default().opposite(), result);
    }

    #[test]
    fn test_calculate_pole_azimuth_and_direction() {
        // Greenwich equator
        let g_eq = Vector3d::new(1.0, 0.0, 0.0);

        // 90 degrees East on the equator
        let e_eq = Vector3d::new(0.0, 1.0, 0.0);

        // 90 degrees West on the equator
        let w_eq = Vector3d::new(0.0, -1.0, 0.0);

        let angle_90 = Angle::from(Degrees(90.0));
        let pole_a = calculate_pole(
            Angle::from(Degrees(0.0)),
            Angle::from(Degrees(0.0)),
            angle_90,
        );
        assert!(are_orthogonal(&g_eq, &pole_a));

        let dir_a = calculate_direction(
            Angle::from(Degrees(0.0)),
            Angle::from(Degrees(0.0)),
            angle_90,
        );
        assert!(are_orthogonal(&g_eq, &dir_a));
        assert!(are_orthogonal(&pole_a, &dir_a));
        assert_eq!(dir_a, direction(&g_eq, &pole_a));

        let north_pole = Vector3d::new(0.0, 0.0, 1.0);
        assert_eq!(north_pole, pole_a);

        let result = g_eq.cross(&e_eq);
        assert_eq!(north_pole, result);

        let result = calculate_azimuth(&g_eq, &pole_a);
        assert_eq!(angle_90, result);

        let pole_b = calculate_pole(
            Angle::from(Degrees(0.0)),
            Angle::from(Degrees(0.0)),
            -angle_90,
        );
        assert!(are_orthogonal(&g_eq, &pole_b));

        let dir_b = calculate_direction(
            Angle::from(Degrees(0.0)),
            Angle::from(Degrees(0.0)),
            -angle_90,
        );
        assert!(are_orthogonal(&g_eq, &dir_b));
        assert!(are_orthogonal(&pole_b, &dir_b));
        assert_eq!(dir_b, direction(&g_eq, &pole_b));

        let south_pole = Vector3d::new(0.0, 0.0, -1.0);
        assert_eq!(south_pole, pole_b);

        let result = g_eq.cross(&w_eq);
        assert_eq!(south_pole, result);

        let result = calculate_azimuth(&g_eq, &pole_b);
        assert_eq!(-angle_90, result);
    }


    #[test]
    fn test_calculate_position() {
        // Greenwich equator
        let g_eq = Vector3d::new(1.0, 0.0, 0.0);

        // 90 degrees East on the equator
        let e_eq = Vector3d::new(0.0, 1.0, 0.0);

        let pole_0 = g_eq.cross(&e_eq);

        let angle_90 = Angle::from(Degrees(90.0));

        let pos_1 = position(&g_eq, &direction(&g_eq, &pole_0), angle_90);
        assert_eq!(e_eq, pos_1);

        let pos_2 = rotate_position(&g_eq, &pole_0, Angle::default(), angle_90);
        assert_eq!(e_eq, pos_2);

        let pos_3 = rotate_position(&g_eq, &pole_0, angle_90, angle_90);
        assert_eq!(pole_0, pos_3);
    }

    #[test]
    fn test_calculate_cross_track_distance_and_square() {
        // Greenwich equator
        let g_eq = Vector3d::new(1.0, 0.0, 0.0);

        // 90 degrees East on the equator
        let e_eq = Vector3d::new(0.0, 1.0, 0.0);

        let pole_0 = g_eq.cross(&e_eq);

        let longitude = Degrees(1.0);

        // Accuracy drops off outside of this range
        for lat in -83..84 {
            let latitude = Degrees(lat as f64);
            let latlong = LatLong::new(latitude, longitude);
            let point = Vector3d::from(&latlong);

            let expected = (lat as f64).to_radians();
            let xtd = cross_track_distance(&pole_0, &point);
            assert!(is_within_tolerance(
                expected,
                xtd.0,
                2.0 * core::f64::EPSILON
            ));

            let expected = great_circle::gc2e_distance(Radians(expected));
            let expected = expected * expected;
            let xtd2 = sq_cross_track_distance(&pole_0, &point);
            assert!(is_within_tolerance(
                expected,
                xtd2,
                4.0 * core::f64::EPSILON
            ));
        }
    }

    #[test]
    fn test_calculate_along_track_distance_and_square() {
        // Greenwich equator
        let g_eq = Vector3d::new(1.0, 0.0, 0.0);

        // 90 degrees East on the equator
        let e_eq = Vector3d::new(0.0, 1.0, 0.0);

        let pole_0 = g_eq.cross(&e_eq);

        let latitude = Degrees(1.0);

        // Accuracy drops off outside of this range
        for lon in 123..124 {
            let longitude = Degrees(lon as f64);
            let latlong = LatLong::new(latitude, longitude);
            let point = Vector3d::from(&latlong);

            let expected = (lon as f64).to_radians();
            let atd = along_track_distance(&g_eq, &pole_0, &point);
            assert!(is_within_tolerance(
                expected,
                atd.0,
                2.0 * core::f64::EPSILON
            ));

            let (atd, xtd) = calculate_atd_and_xtd(&g_eq, &pole_0, &point);
            assert!(is_within_tolerance(
                expected,
                atd.0,
                2.0 * core::f64::EPSILON
            ));
            assert!(is_within_tolerance(
                1_f64.to_radians(),
                xtd.0,
                2.0 * core::f64::EPSILON
            ));

            let expected = great_circle::gc2e_distance(Radians(expected));
            let expected = expected * expected;
            let atd2 = sq_along_track_distance(&g_eq, &pole_0, &point);
            assert!(is_within_tolerance(
                expected,
                atd2,
                2.0 * core::f64::EPSILON
            ));
        }
    }

    #[test]
    fn test_special_cases() {
        // Greenwich equator
        let g_eq = Vector3d::new(1.0, 0.0, 0.0);

        // 90 degrees East on the equator
        let e_eq = Vector3d::new(0.0, 1.0, 0.0);

        let pole_0 = g_eq.cross(&e_eq);

        assert_eq!(0.0, along_track_distance(&g_eq, &pole_0, &g_eq).0);
        assert_eq!(0.0, sq_along_track_distance(&g_eq, &pole_0, &pole_0));

        let (atd, xtd) = calculate_atd_and_xtd(&g_eq, &pole_0, &g_eq);
        assert_eq!(0.0, atd.0);
        assert_eq!(0.0, xtd.0);

        let (atd, xtd) = calculate_atd_and_xtd(&g_eq, &pole_0, &pole_0);
        assert_eq!(0.0, atd.0);
        assert_eq!(core::f64::consts::FRAC_PI_2, xtd.0);

        let (atd, xtd) = calculate_atd_and_xtd(&g_eq, &pole_0, &-pole_0);
        assert_eq!(0.0, atd.0);
        assert_eq!(-core::f64::consts::FRAC_PI_2, xtd.0);

        // Test for 100% code coverage
        let near_north_pole = LatLong::new(Degrees(89.99999), Degrees(0.0));
        let p = Vector3d::from(&near_north_pole);
        let (atd, xtd) = calculate_atd_and_xtd(&g_eq, &pole_0, &p);
        assert_eq!(0.0, atd.0);
        assert!(is_within_tolerance(
            core::f64::consts::FRAC_PI_2,
            xtd.0,
            0.000001
        ));
    }
}
