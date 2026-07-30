// Copyright (c) 2024-2026 Ken Barker

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
//! calculations using `Vector3`s to represent points and great circle poles
//! on a unit sphere.
//!
//! A `Vector3` is a [nalgebra](https://crates.io/crates/nalgebra) `Vector3`.

extern crate nalgebra as na;

use crate::{Vector3, great_circle};
use angle_sc::{Angle, Radians, trig};
use num_traits::{Float, float::FloatConst};

pub mod intersection;

/// The minimum value of the sine of an f64 angle to normalise.
/// Approximately 7.504e-9 seconds
pub const MIN_SIN_MULTIPLE: u32 = 16384;

/// The minimum value of the sine of an f32 angle to normalise.
pub const MIN_SIN_MULTIPLE_F32: u32 = 4096;

/// Convert a latitude and longitude to a point on the unit sphere.
///
/// @pre |lat| <= 90.0 degrees.
/// * `lat` - the latitude.
/// * `lon` - the longitude.
///
/// returns a `Vector3` of the point on the unit sphere.
#[must_use]
pub fn to_point<T: Float>(lat: Angle<T>, lon: Angle<T>) -> Vector3<T> {
    Vector3::<T>::new(
        lat.cos().0 * lon.cos().0,
        lat.cos().0 * lon.sin().0,
        lat.sin().0,
    )
}

/// Calculate the latitude of a point.
///
/// * `a` - the point.
///
/// returns the latitude of the point
#[must_use]
pub fn latitude<T: Float>(a: &Vector3<T>) -> Angle<T> {
    Angle::from_y_x(a[2], a[0].hypot(a[1]))
}

/// Calculate the longitude of a point.
///
/// * `a` - the point.
///
/// returns the longitude of the point
#[must_use]
pub fn longitude<T: Float>(a: &Vector3<T>) -> Angle<T> {
    Angle::from_y_x(a[1], a[0])
}

/// Determine whether a `Vector3` is a unit vector.
///
/// * `a` - the vector.
///
/// returns true if `a` is a unit vector, false otherwise.
#[allow(clippy::missing_panics_doc)]
#[must_use]
pub fn is_unit<T>(a: &Vector3<T>) -> bool
where
    T: Float + na::Scalar + na::ComplexField<RealField = T>,
{
    let twelve = T::from(12).expect("Could not convert constant to Float");
    let min_sq_length = T::one() - twelve * T::epsilon();
    let max_sq_length = T::one() + twelve * T::epsilon();

    (min_sq_length..=max_sq_length).contains(&(a.norm_squared()))
}

/// Normalize a vector to lie on the surface of the unit sphere.
///
/// Note: this function returns an `Option` so uses the British spelling of
/// `normalise` to differentiate it from the standard `normalize` function.
/// * `a` the `Vector3`
/// * `min_sq_value` the minimum square of a vector length to normalize.
///
/// return the nomalized point or None if the vector is too small to normalize.
#[must_use]
pub fn normalise<T>(a: &Vector3<T>, min_sq_value: T) -> Option<Vector3<T>>
where
    T: Float + na::Scalar + na::ComplexField<RealField = T>,
{
    if a.norm_squared() < min_sq_value {
        None
    } else {
        Some(a.normalize())
    }
}

/// Calculate the square of the Euclidean distance between two points.
///
/// Note: points do NOT need to be valid Points.
/// @post for unit vectors: result <= 4
/// * `a`, `b` the points.
///
/// returns the square of the Euclidean distance between the points.
#[must_use]
pub fn sq_distance<T>(a: &Vector3<T>, b: &Vector3<T>) -> T
where
    T: Float + na::Scalar + na::ComplexField<RealField = T>,
{
    (b - a).norm_squared()
}

/// Calculate the shortest (Euclidean) distance between two Points.
///
/// @post for unit vectors: result <= 2
/// * `a`, `b` the points.
///
/// returns the shortest (Euclidean) distance between the points.
#[must_use]
pub fn distance<T>(a: &Vector3<T>, b: &Vector3<T>) -> T
where
    T: Float + na::Scalar + na::ComplexField<RealField = T>,
{
    (b - a).norm()
}

/// Determine whether two `Vector3`s are orthogonal (perpendicular).
///
/// * `a`, `b` the `Vector3`s.
///
/// returns true if a and b are orthogonal, false otherwise.
#[must_use]
pub fn are_orthogonal<T>(a: &Vector3<T>, b: &Vector3<T>) -> bool
where
    T: Float + na::Scalar + na::ComplexField<RealField = T>,
{
    let two_epsilon = T::epsilon() + T::epsilon();
    let max_length = two_epsilon + two_epsilon;

    (-max_length..=max_length).contains(&(a.dot(b)))
}

/// Calculate the relative longitude of point a from point b.
///
/// * `a`, `b` - the points.
///
/// returns the relative longitude of point a from point b,
/// negative if a is West of b, positive otherwise.
#[must_use]
pub fn delta_longitude<T>(a: &Vector3<T>, b: &Vector3<T>) -> Angle<T>
where
    T: Float + na::Scalar + na::ComplexField<RealField = T>,
{
    let a_lon = a.xy();
    let b_lon = b.xy();
    Angle::from_y_x(b_lon.perp(&a_lon), b_lon.dot(&a_lon))
}

/// Determine whether point a is South of point b.
///
/// It calculates and compares the z component of the two points.
/// * `a`, `b` - the points.
///
/// returns true if a is South of b, false otherwise.
#[must_use]
pub fn is_south_of<T: Float>(a: &Vector3<T>, b: &Vector3<T>) -> bool {
    a[2] < b[2]
}
/// Determine whether point a is West of point b.
///
/// It calculates and compares the perp product of the two points.
/// * `a`, `b` - the points.
///
/// returns true if a is West of b, false otherwise.
#[must_use]
pub fn is_west_of<T>(a: &Vector3<T>, b: &Vector3<T>) -> bool
where
    T: Float + na::Scalar + na::ComplexField<RealField = T>,
{
    b.xy().perp(&a.xy()) < T::zero()
}

/// Calculate the right hand pole vector of a Great Circle from an initial
/// position and an azimuth.
///
/// See: <http://www.movable-type.co.uk/scripts/latlong-vectors.html#distance>
/// * `lat` - start point Latitude.
/// * `lon` - start point Longitude.
/// * `azi` - start point azimuth.
///
/// returns the right hand pole vector of the great circle.
#[must_use]
pub fn calculate_pole<T: Float>(lat: Angle<T>, lon: Angle<T>, azi: Angle<T>) -> Vector3<T> {
    let x = trig::UnitNegRange::<T>::clamp(
        lon.sin().0 * azi.cos().0 - lat.sin().0 * lon.cos().0 * azi.sin().0,
    );
    let y = trig::UnitNegRange::<T>::clamp(
        T::zero() - lon.cos().0 * azi.cos().0 - lat.sin().0 * lon.sin().0 * azi.sin().0,
    );
    let z = trig::UnitNegRange::<T>(lat.cos().0 * azi.sin().0);

    Vector3::new(x.0, y.0, z.0)
}

/// Calculate the azimuth at a point on the Great Circle defined by pole.
///
/// * `point` - the point.
/// * `pole` - the right hand pole of the Great Circle.
///
/// returns the azimuth at the point on the great circle.
#[must_use]
pub fn calculate_azimuth<T>(point: &Vector3<T>, pole: &Vector3<T>) -> Angle<T>
where
    T: Float + na::Scalar + na::ComplexField<RealField = T>,
{
    let max_lat = T::one() - (T::epsilon() + T::epsilon());

    let sin_lat: T = point[2];
    // if the point is close to the North or South poles, azimuth is 180 or 0.
    if max_lat <= Float::abs(sin_lat) {
        // azimuth is zero or 180 degrees
        return if sin_lat.is_sign_negative() {
            Angle::default()
        } else {
            Angle::new(trig::UnitNegRange(T::zero()), trig::UnitNegRange(-T::one()))
        };
    }

    Angle::from_y_x(pole[2], pole.xy().perp(&point.xy()))
}

/// Calculate the direction vector along a Great Circle from an initial
/// position and an azimuth.
///
/// See: Panou and Korakitis equations: 30, 31, & 32a
/// <https://arxiv.org/abs/1811.03513>
/// * `lat` - start point Latitude.
/// * `lon` - start point Longitude.
/// * `azi` - start point azimuth.
///
/// returns the direction vector at the point on the great circle.
#[must_use]
pub fn calculate_direction<T: Float>(lat: Angle<T>, lon: Angle<T>, azi: Angle<T>) -> Vector3<T> {
    let x = trig::UnitNegRange::clamp(
        T::zero() - lat.sin().0 * lon.cos().0 * azi.cos().0 - lon.sin().0 * azi.sin().0,
    );
    let y = trig::UnitNegRange::clamp(
        T::zero() - lat.sin().0 * lon.sin().0 * azi.cos().0 + lon.cos().0 * azi.sin().0,
    );
    let z = trig::UnitNegRange(lat.cos().0 * azi.cos().0);

    Vector3::new(x.0, y.0, z.0)
}

/// Calculate the direction vector of a Great Circle arc.
///
/// * `a` - the start point.
/// * `pole` - the pole of a Great Circle.
///
/// returns the direction vector at the point on the great circle.
#[must_use]
pub fn direction<T>(a: &Vector3<T>, pole: &Vector3<T>) -> Vector3<T>
where
    T: Float + na::Scalar + na::ComplexField<RealField = T>,
{
    pole.cross(a)
}

/// Calculate the position of a point along a Great Circle arc.
///
/// * `a` - the start point.
/// * `dir` - the direction vector of a Great Circle at a.
/// * `distance` - the a Great Circle as an Angle.
///
/// returns the position vector at the point on the great circle.
#[must_use]
pub fn position<T>(a: &Vector3<T>, dir: &Vector3<T>, distance: Angle<T>) -> Vector3<T>
where
    T: Float + na::Scalar + na::ComplexField<RealField = T>,
{
    a * distance.cos().0 + dir * distance.sin().0
}

/// Calculate the direction vector of a Great Circle rotated by angle.
///
/// * `dir` - the direction vector of a Great Circle arc.
/// * `pole` - the pole of a Great Circle.
/// * `angle` - the angle to rotate the direction vector by.
///
/// returns the direction vector at the point on the great circle
/// rotated by angle.
#[must_use]
pub fn rotate<T>(dir: &Vector3<T>, pole: &Vector3<T>, angle: Angle<T>) -> Vector3<T>
where
    T: Float + na::Scalar + na::ComplexField<RealField = T>,
{
    position(dir, pole, angle)
}

/// Calculate the position of a point rotated by angle at radius.
///
/// * `a` - the start point.
/// * `pole` - the pole of a Great Circle.
/// * `angle` - the angle to rotate the direction vector by.
/// * `radius` - the radius from the start point.
///
/// returns the position vector at angle and radius from the start point.
#[must_use]
pub fn rotate_position<T>(
    a: &Vector3<T>,
    pole: &Vector3<T>,
    angle: Angle<T>,
    radius: Angle<T>,
) -> Vector3<T>
where
    T: Float + na::Scalar + na::ComplexField<RealField = T>,
{
    position(a, &rotate(&direction(a, pole), pole, angle), radius)
}

/// The sine of the across track distance of a point relative to a Great Circle pole.
///
/// It is simply the dot product of the pole and the point: pole . point
/// * `pole` - the Great Circle pole.
/// * `point` - the point.
///
/// returns the sine of the across track distance of point relative to the pole.
#[must_use]
fn sin_xtd<T>(pole: &Vector3<T>, point: &Vector3<T>) -> trig::UnitNegRange<T>
where
    T: Float + na::Scalar + na::ComplexField<RealField = T>,
{
    trig::UnitNegRange::clamp(pole.dot(point))
}

/// Determine whether point is right of a Great Circle pole.
///
/// It compares the dot product of the pole and point.
/// * `pole` - the Great Circle pole.
/// * `point` - the point.
///
/// returns true if the point is right of the pole,
/// false if on or to the left of the Great Circle.
#[must_use]
pub fn is_right_of<T>(pole: &Vector3<T>, point: &Vector3<T>) -> bool
where
    T: Float + na::Scalar + na::ComplexField<RealField = T>,
{
    pole.dot(point) < T::zero()
}

/// The across track distance of a point relative to a Great Circle pole.
///
/// * `pole` - the Great Circle pole.
/// * `point` - the point.
///
/// returns the across track distance of point relative to pole, in `Radians`.
#[must_use]
pub fn cross_track_distance<T>(pole: &Vector3<T>, point: &Vector3<T>) -> Radians<T>
where
    T: Float + FloatConst + na::Scalar + na::ComplexField<RealField = T>,
{
    let sin_d = sin_xtd(pole, point);
    if Float::abs(sin_d.0) < T::epsilon() {
        Radians(T::zero())
    } else {
        Radians(Float::asin(sin_d.0))
    }
}

/// The square of the Euclidean cross track distance of a point relative to a
/// Great Circle pole.
///
/// * `pole` - the Great Circle pole.
/// * `point` - the point.
///
/// returns the square of the euclidean distance of point relative to pole.
#[must_use]
pub fn sq_cross_track_distance<T>(pole: &Vector3<T>, point: &Vector3<T>) -> T
where
    T: Float + na::Scalar + na::ComplexField<RealField = T>,
{
    let two = T::one() + T::one();
    let sin_d = sin_xtd(pole, point);
    if Float::abs(sin_d.0) < T::epsilon() {
        T::zero()
    } else {
        two * (T::one() - trig::swap_sin_cos(sin_d).0)
    }
}

/// Calculate the closest point on a plane to the given point.
///
/// See: [Closest Point on Plane](https://gdbooks.gitbooks.io/3dcollisions/content/Chapter1/closest_point_on_plane.html)
/// * `pole` - the Great Circle pole (aka normal) of the plane.
/// * `point` - the point.
///
/// returns the closest point on a plane to the given point.
#[must_use]
fn calculate_point_on_plane<T>(pole: &Vector3<T>, point: &Vector3<T>) -> Vector3<T>
where
    T: Float + na::Scalar + na::ComplexField<RealField = T>,
{
    let t = sin_xtd(pole, point);
    point - pole * t.0
}

/// The sine of the along track distance of a point along a Great Circle arc.
///
/// It is the triple product of the pole, a and the point:
/// (pole X a) . point = pole . (a X point)
/// * `a` - the start point of the Great Circle arc.
/// * `pole` - the pole of the Great Circle arc.
/// * `point` - the point.
///
/// returns the sine of the along track distance of point relative to the start
/// of a great circle arc.
#[must_use]
fn sin_atd<T>(a: &Vector3<T>, pole: &Vector3<T>, point: &Vector3<T>) -> trig::UnitNegRange<T>
where
    T: Float + na::Scalar + na::ComplexField<RealField = T>,
{
    trig::UnitNegRange::clamp(pole.cross(a).dot(point))
}

/// Calculate the relative distance of two points on a Great Circle arc.
///
/// @pre both points must be on the Great Circle defined by `pole`.
/// * `a` - the start point of the Great Circle arc.
/// * `pole` - the pole of the Great Circle arc.
/// * `point` - a point in the Great Circle.
///
/// returns the Great Circle along track distance in `Radians`.
#[must_use]
pub fn calculate_great_circle_atd<T>(
    a: &Vector3<T>,
    pole: &Vector3<T>,
    point: &Vector3<T>,
) -> Radians<T>
where
    T: Float + FloatConst + na::Scalar + na::ComplexField<RealField = T>,
{
    let min_distance = T::epsilon() + T::epsilon();
    let min_sq_distance = min_distance * min_distance;

    let sq_atd = sq_distance(a, point);
    if sq_atd < min_sq_distance {
        Radians(T::zero())
    } else {
        Radians(
            great_circle::e2gc_distance(Float::sqrt(sq_atd))
                .0
                .copysign(sin_atd(a, pole, point).0),
        )
    }
}

/// The Great Circle distance of a point along the arc relative to a,
/// (+ve) ahead of a, (-ve) behind a.
///
/// * `a` - the start point of the Great Circle arc.
/// * `pole` - the pole of the Great Circle arc.
/// * `point` - the point.
///
/// returns the along track distance of point relative to the start of a great circle arc.
#[allow(clippy::missing_panics_doc)]
#[must_use]
pub fn along_track_distance<T>(a: &Vector3<T>, pole: &Vector3<T>, point: &Vector3<T>) -> Radians<T>
where
    T: Float + FloatConst + na::Scalar + na::ComplexField<RealField = T>,
    f64: From<T>,
{
    let min_angle_multiple = if f64::from(T::epsilon()) < f64::epsilon() {
        MIN_SIN_MULTIPLE_F32
    } else {
        MIN_SIN_MULTIPLE
    };
    let min_angle_multiple =
        T::from(min_angle_multiple).expect("Could not convert constant to Float");
    let min_sin_angle = min_angle_multiple * T::epsilon();
    let min_sq_norm = min_sin_angle * min_sin_angle;

    let plane_point = calculate_point_on_plane(pole, point);
    normalise(&plane_point, min_sq_norm).map_or_else(
        || Radians(T::zero()), // point is too close to a pole
        |c| calculate_great_circle_atd(a, pole, &c),
    )
}

/// Calculate the square of the Euclidean along track distance of a point
/// from the start of an Arc.
///
/// It is calculated using the closest point on the plane to the point.
/// * `a` - the start point of the Great Circle arc.
/// * `pole` - the pole of the Great Circle arc.
/// * `point` - the point.
///
/// returns the square of the Euclidean along track distance
#[allow(clippy::missing_panics_doc)]
#[must_use]
pub fn sq_along_track_distance<T>(a: &Vector3<T>, pole: &Vector3<T>, point: &Vector3<T>) -> T
where
    T: Float + na::Scalar + na::ComplexField<RealField = T>,
    f64: From<T>,
{
    let min_distance = T::epsilon() + T::epsilon();
    let min_sq_distance = min_distance * min_distance;

    let min_angle_multiple = if f64::from(T::epsilon()) < f64::epsilon() {
        MIN_SIN_MULTIPLE_F32
    } else {
        MIN_SIN_MULTIPLE
    };
    let min_angle_multiple =
        T::from(min_angle_multiple).expect("Could not convert constant to Float");
    let min_sin_angle = min_angle_multiple * T::epsilon();
    let min_sq_norm = min_sin_angle * min_sin_angle;

    let plane_point = calculate_point_on_plane(pole, point);
    normalise(&plane_point, min_sq_norm).map_or_else(
        || T::zero(), // point is too close to a pole
        |c| {
            let sq_d = sq_distance(a, &(c));
            if sq_d < min_sq_distance {
                T::zero()
            } else {
                sq_d
            }
        },
    )
}

/// Calculate Great Circle along and across track distances.
///
/// * `a` - the start point of the Great Circle arc.
/// * `pole` - the pole of the Great Circle arc.
/// * `p` - the point.
///
/// returns the along and across track distances of point relative to the
/// start of a great circle arc.
#[allow(clippy::missing_panics_doc)]
#[allow(clippy::similar_names)]
#[must_use]
pub fn calculate_atd_and_xtd<T>(
    a: &Vector3<T>,
    pole: &Vector3<T>,
    p: &Vector3<T>,
) -> (Radians<T>, Radians<T>)
where
    T: Float + FloatConst + na::Scalar + na::ComplexField<RealField = T>,
    f64: From<T>,
{
    let min_distance = T::epsilon() + T::epsilon();
    let min_sq_distance = min_distance * min_distance;

    let min_angle_multiple = if f64::from(T::epsilon()) < f64::epsilon() {
        MIN_SIN_MULTIPLE_F32
    } else {
        MIN_SIN_MULTIPLE
    };
    let min_angle_multiple =
        T::from(min_angle_multiple).expect("Could not convert constant to Float");
    let min_sin_angle = min_angle_multiple * T::epsilon();
    let min_sq_norm = min_sin_angle * min_sin_angle;

    let mut atd = Radians(T::zero());
    let mut xtd = Radians(T::zero());

    let sq_d = sq_distance(a, p);
    if sq_d >= min_sq_distance {
        // point is not close to a
        let sin_xtd = sin_xtd(pole, p).0;
        if Float::abs(sin_xtd) >= T::epsilon() {
            xtd = Radians(Float::asin(sin_xtd));
        }

        // the closest point on the plane of the pole to the point
        let plane_point = p - pole * sin_xtd;
        atd = normalise(&plane_point, min_sq_norm).map_or_else(
            || Radians(T::zero()), // point is too close to a pole
            |c| calculate_great_circle_atd(a, pole, &c),
        );
    }

    (atd, xtd)
}

/// Normalise a centroid on coincident great circles.
///
/// Note: it handles the case where the centroid is too small to normalise.
///
/// * `centroid` - the centroid to be normalised.
/// * `point` - a mid-point.
/// * `pole` - the pole of the Great Circle arc.
///
/// returns the normalise centroid.
#[allow(clippy::missing_panics_doc)]
#[must_use]
pub fn normalise_centroid<T>(
    centroid: &Vector3<T>,
    point: &Vector3<T>,
    pole: &Vector3<T>,
) -> Vector3<T>
where
    T: Float + na::Scalar + na::ComplexField<RealField = T>,
    f64: From<T>,
{
    let min_angle_multiple = if f64::from(T::epsilon()) < f64::epsilon() {
        MIN_SIN_MULTIPLE_F32
    } else {
        MIN_SIN_MULTIPLE
    };
    let min_angle_multiple =
        T::from(min_angle_multiple).expect("Could not convert constant to Float");
    let min_sin_angle = min_angle_multiple * T::epsilon();
    let min_sq_norm = min_sin_angle * min_sin_angle;

    normalise(centroid, min_sq_norm).unwrap_or_else(|| {
        // centroid is half way between points

        // calculate a point on the coincident great circle
        // half way between points, closer to the start of the arc
        position(
            point,
            &direction(point, pole),
            Angle::default().quarter_turn_ccw(),
        )
    })
}

#[cfg(test)]
mod tests {
    use core::f64;

    use super::*;
    use crate::LatLong;
    use angle_sc::{Degrees, Radians, is_within_tolerance};

    pub const MIN_SIN_ANGLE: f64 = (MIN_SIN_MULTIPLE as f64) * f64::EPSILON;
    pub const MIN_SQ_NORM: f64 = MIN_SIN_ANGLE * MIN_SIN_ANGLE;

    #[test]
    fn test_normalise() {
        let zero = Vector3::new(0.0, 0.0, 0.0);
        assert!(normalise(&zero, MIN_SQ_NORM).is_none());

        // Greenwich equator
        let g_eq = Vector3::new(1.0, 0.0, 0.0);
        assert!(normalise(&g_eq, MIN_SQ_NORM).is_some());

        // A vector just too small to normalize
        let too_small = Vector3::new(16383.0 * f64::EPSILON, 0.0, 0.0);
        assert!(normalise(&too_small, MIN_SQ_NORM).is_none());

        assert_eq!(2.0844083160439303e-10, MIN_SIN_ANGLE.to_degrees().asin());

        // A vector just large enough to normalize
        let small = Vector3::new(MIN_SIN_ANGLE, 0.0, 0.0);
        let result = normalise(&small, MIN_SQ_NORM);
        assert!(result.is_some());

        assert!(is_unit(&result.unwrap()));
        assert_eq!(result.unwrap(), g_eq);
    }

    #[test]
    fn test_point_lat_longs() {
        // Test South pole
        let lat_lon_south = LatLong::new(Degrees(-90.0), Degrees(180.0));
        let point_south = Vector3::from(&lat_lon_south);
        assert!(is_unit(&point_south));
        assert_eq!(Vector3::new(0.0, 0.0, -1.0), point_south);

        assert_eq!(Degrees(-90.0), Degrees::from(latitude(&point_south)));
        assert_eq!(Degrees(0.0), Degrees::from(longitude(&point_south)));

        let result = LatLong::from(&point_south);
        assert_eq!(-90.0, result.lat().0);
        // Note: longitude is now zero, since the poles do not have a Longitude
        assert_eq!(0.0, result.lon().0);

        // Test Greenwich equator
        let lat_lon_0_0 = LatLong::new(Degrees(0.0), Degrees(0.0));
        let point_0 = Vector3::from(&lat_lon_0_0);
        assert!(is_unit(&point_0));
        assert_eq!(Vector3::new(1.0, 0.0, 0.0), point_0);
        assert_eq!(lat_lon_0_0, LatLong::from(&point_0));

        // Test antimeridian equator
        let lat_lon_0_180 = LatLong::new(Degrees(0.0), Degrees(180.0));
        let point_1 = Vector3::from(&lat_lon_0_180);
        assert!(is_unit(&point_1));
        assert_eq!(Vector3::new(-1.0, 0.0, 0.0), point_1);
        assert_eq!(false, is_west_of(&point_0, &point_1));
        assert_eq!(
            Radians(core::f64::consts::PI),
            Radians::from(delta_longitude(&point_0, &point_1)).abs()
        );

        let lat_lon_0_m180 = LatLong::new(Degrees(0.0), Degrees(-180.0));
        let point_2 = Vector3::from(&lat_lon_0_m180);
        assert!(is_unit(&point_2));
        assert_eq!(Vector3::new(-1.0, 0.0, 0.0), point_2);
        // Converts back to +ve longitude
        assert_eq!(lat_lon_0_180, LatLong::from(&point_2));

        assert_eq!(false, is_west_of(&point_0, &point_2));
        assert_eq!(
            -core::f64::consts::PI,
            Radians::from(delta_longitude(&point_0, &point_2)).0
        );

        let lat_lon_0_r3 = LatLong::new(Degrees(0.0), Degrees(3.0_f64.to_degrees()));
        let point_3 = Vector3::from(&lat_lon_0_r3);
        assert!(is_unit(&point_3));
        let result = LatLong::from(&point_3);
        assert_eq!(0.0, result.lat().0);
        assert_eq!(
            3.0_f64,
            Radians::from(delta_longitude(&point_3, &point_0)).0
        );
        assert_eq!(3.0_f64.to_degrees(), result.lon().0);
        assert!(is_west_of(&point_0, &point_3));
        assert_eq!(-3.0, Radians::from(delta_longitude(&point_0, &point_3)).0);

        assert_eq!(false, is_west_of(&point_1, &point_3));
        assert!(is_within_tolerance(
            core::f64::consts::PI - 3.0,
            Radians::from(delta_longitude(&point_1, &point_3)).0,
            f64::EPSILON
        ));

        let lat_lon_0_mr3 = LatLong::new(Degrees(0.0), Degrees(-3.0_f64.to_degrees()));
        let point_4 = Vector3::from(&lat_lon_0_mr3);
        assert!(is_unit(&point_4));
        assert_eq!(3.0, Radians::from(delta_longitude(&point_0, &point_4)).0);

        let result = LatLong::from(&point_4);
        assert_eq!(0.0, result.lat().0);
        assert_eq!(-3.0_f64.to_degrees(), result.lon().0);
        assert!(is_west_of(&point_1, &point_4));
        assert!(is_within_tolerance(
            3.0 - core::f64::consts::PI,
            Radians::from(delta_longitude(&point_1, &point_4)).0,
            f64::EPSILON
        ));
    }

    #[test]
    fn test_point_distance() {
        let lat_lon_south = LatLong::new(Degrees(-90.0), Degrees(0.0));
        let south_pole = Vector3::from(&lat_lon_south);

        let lat_lon_north = LatLong::new(Degrees(90.0), Degrees(0.0));
        let north_pole = Vector3::from(&lat_lon_north);

        assert_eq!(0.0, sq_distance(&south_pole, &south_pole));
        assert_eq!(0.0, sq_distance(&north_pole, &north_pole));
        assert_eq!(4.0, sq_distance(&south_pole, &north_pole));

        assert_eq!(0.0, distance(&south_pole, &south_pole));
        assert_eq!(0.0, distance(&north_pole, &north_pole));
        assert_eq!(2.0, distance(&south_pole, &north_pole));

        // Greenwich equator
        let g_eq = Vector3::new(1.0, 0.0, 0.0);

        // Test IDL equator
        let idl_eq = Vector3::new(-1.0, 0.0, 0.0);

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
        let g_eq = Vector3::new(1.0, 0.0, 0.0);
        let south_pole = Vector3::new(0.0, 0.0, -1.0);
        let result = calculate_azimuth(&south_pole, &g_eq);
        assert_eq!(Angle::default(), result);

        let north_pole = Vector3::new(0.0, 0.0, 1.0);
        let result = calculate_azimuth(&north_pole, &g_eq);
        assert_eq!(Angle::default().opposite(), result);
    }

    #[test]
    fn test_calculate_pole_azimuth_and_direction() {
        // Greenwich equator
        let g_eq = Vector3::new(1.0, 0.0, 0.0);

        // 90 degrees East on the equator
        let e_eq = Vector3::new(0.0, 1.0, 0.0);

        // 90 degrees West on the equator
        let w_eq = Vector3::new(0.0, -1.0, 0.0);

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

        let north_pole = Vector3::new(0.0, 0.0, 1.0);
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

        let south_pole = Vector3::new(0.0, 0.0, -1.0);
        assert_eq!(south_pole, pole_b);

        let result = g_eq.cross(&w_eq);
        assert_eq!(south_pole, result);

        let result = calculate_azimuth(&g_eq, &pole_b);
        assert_eq!(-angle_90, result);
    }

    #[test]
    fn test_calculate_position() {
        // Greenwich equator
        let g_eq = Vector3::new(1.0, 0.0, 0.0);

        // 90 degrees East on the equator
        let e_eq = Vector3::new(0.0, 1.0, 0.0);

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
        let g_eq = Vector3::new(1.0, 0.0, 0.0);

        // 90 degrees East on the equator
        let e_eq = Vector3::new(0.0, 1.0, 0.0);

        let pole_0 = g_eq.cross(&e_eq);

        let longitude = Degrees(1.0);

        for lat in -89..90 {
            let latitude = Degrees(f64::from(lat));
            let latlong = LatLong::new(latitude, longitude);
            let point = Vector3::from(&latlong);

            assert_eq!(lat < 0, is_south_of(&point, &g_eq));
            assert_eq!(lat >= 0, !is_south_of(&point, &e_eq));
            assert_eq!(lat < 0, is_right_of(&pole_0, &point));

            let expected = (f64::from(lat)).to_radians();
            let xtd = cross_track_distance(&pole_0, &point);
            // Accuracy reduces outside of this range
            let tolerance = if (-83..84).contains(&lat) {
                2.0 * f64::EPSILON
            } else {
                32.0 * f64::EPSILON
            };
            assert!(is_within_tolerance(expected, xtd.0, tolerance));

            let expected = great_circle::gc2e_distance(Radians(expected));
            let expected = expected * expected;
            let xtd2 = sq_cross_track_distance(&pole_0, &point);
            // Accuracy reduces outside of this range
            let tolerance = if (-83..84).contains(&lat) {
                4.0 * f64::EPSILON
            } else {
                64.0 * f64::EPSILON
            };
            assert!(is_within_tolerance(expected, xtd2, tolerance));
        }
    }

    #[test]
    fn test_calculate_along_track_distance_and_square() {
        // Greenwich equator
        let g_eq = Vector3::new(1.0, 0.0, 0.0);

        // 90 degrees East on the equator
        let e_eq = Vector3::new(0.0, 1.0, 0.0);

        let pole_0 = g_eq.cross(&e_eq);

        // North of Equator
        let latitude = Degrees(1.0);

        for lon in -179..180 {
            let longitude = Degrees(f64::from(lon));
            let latlong = LatLong::new(latitude, longitude);
            let point = Vector3::from(&latlong);

            let expected = (f64::from(lon)).to_radians();
            let atd = along_track_distance(&g_eq, &pole_0, &point);
            // Accuracy reduces outside of this range
            let tolerance = if (-153..154).contains(&lon) {
                4.0 * f64::EPSILON
            } else {
                32.0 * f64::EPSILON
            };
            assert!(is_within_tolerance(expected, atd.0, tolerance));

            let (atd, xtd) = calculate_atd_and_xtd(&g_eq, &pole_0, &point);
            assert!(is_within_tolerance(expected, atd.0, tolerance));
            assert!(is_within_tolerance(1_f64.to_radians(), xtd.0, f64::EPSILON));

            let expected = great_circle::gc2e_distance(Radians(expected));
            let expected = expected * expected;
            let atd2 = sq_along_track_distance(&g_eq, &pole_0, &point);
            // Accuracy reduces outside of this range
            let tolerance = if (-86..87).contains(&lon) {
                2.0 * f64::EPSILON
            } else {
                32.0 * f64::EPSILON
            };
            assert!(is_within_tolerance(expected, atd2, tolerance));
        }
    }

    #[test]
    fn test_special_cases() {
        // Greenwich equator
        let g_eq = Vector3::new(1.0, 0.0, 0.0);

        // 90 degrees East on the equator
        let e_eq = Vector3::new(0.0, 1.0, 0.0);

        let pole_0 = g_eq.cross(&e_eq);

        // points are at the poles, so atc and sq_atd are zero
        assert_eq!(0.0, along_track_distance(&g_eq, &pole_0, &pole_0).0);
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
        let p = Vector3::from(&near_north_pole);
        let (atd, xtd) = calculate_atd_and_xtd(&g_eq, &pole_0, &p);
        assert_eq!(0.0, atd.0);
        assert!(is_within_tolerance(
            core::f64::consts::FRAC_PI_2,
            xtd.0,
            0.000001
        ));
    }

    #[test]
    fn test_normalise_centroid() {
        let point_0 = Vector3::new(0.0, 0.0, 0.0);
        let point_1 = Vector3::new(1.0, 0.0, 0.0);
        let point_m1 = -point_1;
        let pole_1 = Vector3::new(0.0, 0.0, 1.0);

        // normalised centroid from point_1
        let result = normalise_centroid(&point_0, &point_1, &pole_1);
        assert_eq!(Vector3::new(0.0, -1.0, 0.0), result);

        // normalised centroid from point_1 antoipodal point
        let result = normalise_centroid(&point_0, &point_m1, &pole_1);
        assert_eq!(Vector3::new(0.0, 1.0, 0.0), result);

        // normalised centroid from point_1 centroid
        let point_2 = point_1 + point_1;
        let result = normalise_centroid(&point_2, &point_1, &pole_1);
        assert_eq!(point_1, result);
    }
}
