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

//! The `intersection` module contains functions for calculating great-circle
//! intersections using vectors.
//!
//! A pair of great circles intersect at two points unless they are coincident.\
//! For example, points `u` and `v` in *Figure1*.
//!
//! ![great circle path](https://upload.wikimedia.org/wikipedia/commons/thumb/c/cb/Illustration_of_great-circle_distance.svg/220px-Illustration_of_great-circle_distance.svg.png)\
//! *Figure 1 A pair of intersecting great circles*
//!
//! A great circle intersection point can simply be calculated by normalizing
//! the [cross product](https://en.wikipedia.org/wiki/Cross_product) of their
//! pole vectors.\
//! If the resulting vector is too small to normalize, then the great circles
//! are coincident, in which case they effectively *intersect* everywhere.
//!
//! If a pair of `Arc`s are on coincident great circles,
//! `calculate_coincident_arc_distances` calculates the distances between
//! `Arc` ends, zero if the `Arc`s overlap.
//!
//! Otherwise `use_antipodal_point` determines which intersection point
//! is closer to the [centroid](https://en.wikipedia.org/wiki/Centroid)
//! of the `Arc`s midpoints.

use super::{
    MIN_SQ_NORM, Vector3d, calculate_great_circle_atd, normalise, normalise_centroid, sq_distance,
};
use angle_sc::{Angle, Radians};

/// Calculate an intersection point between the poles of two Great Circles.
/// See: <http://www.movable-type.co.uk/scripts/latlong-vectors.html#intersection>
/// * `pole_0`, `pole_1` the poles.
/// * `min_sq_value` the minimum square of a vector length to normalize.
///
/// return an intersection point or None if the poles represent coincident Great Circles.
#[must_use]
pub fn calculate_intersection(
    pole_0: &Vector3d,
    pole_1: &Vector3d,
    min_sq_value: f64,
) -> Option<Vector3d> {
    normalise(&pole_0.cross(pole_1), min_sq_value)
}

/// Determine whether the antipodal point is closer to the centroid of the
/// `Arc`s.
///
/// * `point` a great-circle intersection point.
/// * `centroid` the centroid of the `Arc`s mid points.
///
/// returns true if the antipodal intersection is closer to the `centroid`
/// of the `Arc`s otherwise returns false.
#[must_use]
pub fn use_antipodal_point(point: &Vector3d, centroid: &Vector3d) -> bool {
    sq_distance(centroid, &(-*point)) < sq_distance(centroid, point)
}

/// Return the closer intersection point to the centroid of the `Arc`s.
///
/// * `point` a great-circle intersection point.
/// * `centroid` the centroid of the `Arc`s mid points.
///
/// returns the antipodal point if it is closer to the `centroid`,
/// otherwise returns the point.
#[must_use]
pub fn closest_intersection_point(point: &Vector3d, centroid: &Vector3d) -> Vector3d {
    if use_antipodal_point(point, centroid) {
        -*point
    } else {
        *point
    }
}

/// Determine the reference point of a pair of arcs.
/// I.e. the closest intersection point if they intersect or the
/// centroid normalized to lie on the unit sphere if they don't.
///
/// * `mid_point_0`, `mid_point_1` the mid points of the `Arc`s.
/// * `pole_0`, `pole_1` the poles of the `Arc` great circles.
///
/// returns the closest intersection point or normalized centroid and the
/// sine of the angle between the arcs, zero if the arcs are coincident.
/// And the absolute relative angle at the intersection point or centroid.
#[must_use]
pub fn calculate_reference_point_and_angle(
    mid_point_0: &Vector3d,
    pole_0: &Vector3d,
    mid_point_1: &Vector3d,
    pole_1: &Vector3d,
) -> (Vector3d, Angle) {
    let centroid = mid_point_0 + mid_point_1;
    // calculate the intersection point between the great circles
    let point = pole_0.cross(pole_1);
    normalise(&point, MIN_SQ_NORM).map_or_else(
        || {
            // the great circles are coincident

            let c = normalise_centroid(&centroid, mid_point_0, pole_0);

            // determine whether the great circles oppose each other
            let angle = if pole_0.dot(pole_1).is_sign_negative() {
                Angle::default().opposite()
            } else {
                Angle::default()
            };

            (c, angle)
        },
        |p| {
            // the great circles intersect

            let x = closest_intersection_point(&p, &centroid);
            (x, Angle::from_y_x(point.norm(), pole_0.dot(pole_1)))
        },
    )
}

/// Calculate signed great circle distances from two arc mid points to their
/// closest intersection point or normalized centroid if the arcs are on coincident
/// great circles.
///
/// * `mid_point_0`, `mid_point_1` the mid points of the arcs.
/// * `pole_0`, `pole_1` the poles of the arc great circles.
///
/// returns the signed great circle distances of the closest intersection
/// point or centroid  from the arc mid points in `Radians`,
/// and the relative angle between the arc great circles.
#[must_use]
pub fn calculate_arc_reference_distances_and_angle(
    mid_point_0: &Vector3d,
    pole_0: &Vector3d,
    mid_point_1: &Vector3d,
    pole_1: &Vector3d,
) -> (Radians, Radians, Angle) {
    let (point, angle) =
        calculate_reference_point_and_angle(mid_point_0, pole_0, mid_point_1, pole_1);
    let distance_0 = calculate_great_circle_atd(mid_point_0, pole_0, &point);
    let distance_1 = calculate_great_circle_atd(mid_point_1, pole_1, &point);

    (distance_0, distance_1, angle)
}

#[cfg(test)]
mod tests {
    use std::f64;

    use super::*;
    use crate::{LatLong, vector};
    use angle_sc::{Angle, Degrees, is_within_tolerance};

    #[test]
    fn test_calculate_intersection() {
        let lat_lon_south = LatLong::new(Degrees(-90.0), Degrees(0.0));
        let south_pole = Vector3d::from(&lat_lon_south);

        let lat_lon_north = LatLong::new(Degrees(90.0), Degrees(0.0));
        let north_pole = Vector3d::from(&lat_lon_north);

        let lat_lon_idl = LatLong::new(Degrees(0.0), Degrees(180.0));
        let idl = Vector3d::from(&lat_lon_idl);

        let equator_intersection = calculate_intersection(&south_pole, &north_pole, MIN_SQ_NORM);
        assert!(equator_intersection.is_none());

        let gc_intersection1 = calculate_intersection(&idl, &north_pole, MIN_SQ_NORM).unwrap();
        let gc_intersection2 = calculate_intersection(&idl, &south_pole, MIN_SQ_NORM).unwrap();

        assert_eq!(gc_intersection1, -gc_intersection2);
    }

    #[test]
    fn test_calculate_arc_reference_distances_and_angle_coincident_great_circles() {
        let point_1 = Vector3d::new(1.0, 0.0, 0.0);
        let pole_1 = Vector3d::new(0.0, 0.0, 1.0);

        // same mid points and great circles
        let result =
            calculate_arc_reference_distances_and_angle(&point_1, &pole_1, &point_1, &pole_1);
        assert_eq!(Radians(0.0), result.0);
        assert_eq!(Radians(0.0), result.1);
        assert_eq!(Degrees(0.0), Degrees::from(result.2));

        // opposite mid points and same great circles
        let point_m1 = -point_1;
        let result =
            calculate_arc_reference_distances_and_angle(&point_1, &pole_1, &point_m1, &pole_1);
        assert!(is_within_tolerance(
            -f64::consts::FRAC_PI_2,
            result.0.0,
            f64::EPSILON
        ));
        assert!(is_within_tolerance(
            f64::consts::FRAC_PI_2,
            result.1.0,
            f64::EPSILON
        ));
        assert_eq!(Degrees(0.0), Degrees::from(result.2));

        // opposite mid points and great circles
        let pole_m1 = -pole_1;
        let result =
            calculate_arc_reference_distances_and_angle(&point_1, &pole_1, &point_m1, &pole_m1);
        assert!(is_within_tolerance(
            -f64::consts::FRAC_PI_2,
            result.0.0,
            f64::EPSILON
        ));
        assert!(is_within_tolerance(
            -f64::consts::FRAC_PI_2,
            result.1.0,
            f64::EPSILON
        ));
        assert_eq!(Degrees(180.0), Degrees::from(result.2));
    }

    #[test]
    fn test_calculate_arc_reference_distances_and_angle_intersecting_great_circles() {
        let point_1 = Vector3d::new(1.0, 0.0, 0.0);
        let pole_1 = Vector3d::new(0.0, 0.0, 1.0);
        let pole_2 = Vector3d::new(0.0, 1.0, 0.0);

        // intersection, same mid points
        let result =
            calculate_arc_reference_distances_and_angle(&point_1, &pole_1, &point_1, &pole_2);
        assert_eq!(Radians(0.0), result.0);
        assert_eq!(Radians(0.0), result.1);
        assert_eq!(Degrees(90.0), Degrees::from(result.2));

        // intersection, same mid points, acute angle
        let pole_3 = (pole_1 + pole_2).normalize();
        let result =
            calculate_arc_reference_distances_and_angle(&point_1, &pole_1, &point_1, &pole_3);
        assert_eq!(Radians(0.0), result.0);
        assert_eq!(Radians(0.0), result.1);
        assert_eq!(Degrees(45.0), Degrees::from(result.2));

        // intersection, same mid points, obtuse angle
        let pole_m3 = -pole_3;
        let result =
            calculate_arc_reference_distances_and_angle(&point_1, &pole_1, &point_1, &pole_m3);
        assert_eq!(Radians(0.0), result.0);
        assert_eq!(Radians(0.0), result.1);
        assert_eq!(Degrees(135.0), Degrees::from(result.2));

        // intersection, different mid points, acute angle
        let point_2 = vector::position(
            &point_1,
            &vector::direction(&point_1, &pole_3),
            Angle::default().quarter_turn_cw(),
        );
        let result =
            calculate_arc_reference_distances_and_angle(&point_1, &pole_1, &point_2, &pole_3);
        assert_eq!(Radians(0.0), result.0);
        assert!(is_within_tolerance(
            -f64::consts::FRAC_PI_2,
            result.1.0,
            f64::EPSILON
        ));
        assert_eq!(Degrees(45.0), Degrees::from(result.2));

        // intersection, different mid points, obtuse angle
        let result =
            calculate_arc_reference_distances_and_angle(&point_1, &pole_1, &point_2, &pole_m3);
        assert_eq!(Radians(0.0), result.0);
        assert!(is_within_tolerance(
            f64::consts::FRAC_PI_2,
            result.1.0,
            f64::EPSILON
        ));
        assert_eq!(Degrees(135.0), Degrees::from(result.2));
    }
}
