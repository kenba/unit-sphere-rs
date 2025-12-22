// Copyright (c) 2024-2025 Ken Barker

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
//! `calculate_intersection_distances` then calculates great-circle distances
//! along the `Arc`s to the intersection point.

use super::{
    MIN_SQ_DISTANCE, MIN_SQ_NORM, Vector3d, calculate_great_circle_atd, normalise,
    normalise_centroid, sq_distance,
};
use angle_sc::{Angle, Radians, max};

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

/// Calculate the great circle distances to an intersection point from
/// points on a pair of great circle arcs, on different great circles.
///
/// * `a_0`, `a_1` the points on the great circle arcs
/// * `pole_0`, `pole_1` the poles of the great circle arcs
/// * `c` the intersection point
///
/// returns a pair of great circle distances along the arcs to the
/// intersection point in `Radians`.
#[must_use]
pub fn calculate_intersection_distances(
    a_0: &Vector3d,
    pole_0: &Vector3d,
    a_1: &Vector3d,
    pole_1: &Vector3d,
    c: &Vector3d,
) -> (Radians, Radians) {
    (
        calculate_great_circle_atd(a_0, pole_0, c),
        calculate_great_circle_atd(a_1, pole_1, c),
    )
}

/// Whether an along track distance is within an `Arc` length including tolerance.
///
/// * `distance` - the along track distance from the start of the `Arc`.
/// * `length` the length of the `Arc`.
/// * `tolerance` the distance tolerance.
///
/// return true if the along track distance is within the length including tolerance,
/// false otherwise.
#[must_use]
pub fn is_alongside(distance: Radians, length: Radians, tolerance: Radians) -> bool {
    (-tolerance <= distance) && (distance <= length + tolerance)
}

/// Whether an intersection point is within an `Arc`.
///
/// * `distance` - the along track distance to the point from the start of the `Arc`.
/// * `length` the length of the `Arc`.
///
/// return true if the intersection point is within the `Arc`, false otherwise.
#[must_use]
pub fn is_within(distance: f64, length: f64) -> bool {
    (-f64::EPSILON <= distance) && (distance <= length + (f64::EPSILON * (1.0 + length)))
}

/// Calculate the great-circle distances along a pair of `Arc`s on coincident
/// Great Circles to their closest (reference) points.
///
/// * `gc_d` the great-circle distance between the arc start points.
/// * `reciprocal` whether the arcs are in reciprocal directions.
/// * `arc_length_0`, `arc_length_1` the `Arc` lengths in `Radians`.
///
/// returns the distances along the first `Arc` and second `Arc` to their closest
/// (reference) points in `Radians`.
#[must_use]
pub fn calculate_coincident_arc_distances(
    gc_d: Radians,
    reciprocal: bool,
    arc_length_0: Radians,
    arc_length_1: Radians,
) -> (Radians, Radians) {
    if reciprocal {
        // if the arcs intersect
        if is_alongside(
            gc_d,
            max(arc_length_0, arc_length_1),
            Radians(4.0 * f64::EPSILON),
        ) {
            if gc_d <= arc_length_1 {
                // The start of the first `Arc` is within the second `Arc`
                (Radians(0.0), gc_d.clamp(arc_length_1))
            } else {
                // The start of the second `Arc` is within the first `Arc`
                (gc_d.clamp(arc_length_0), Radians(0.0))
            }
        } else {
            let abs_d = gc_d.abs();

            // The distance between the `Arc` b ends
            let b_d = abs_d.0 - arc_length_0.0 - arc_length_1.0;
            // The distance between the `Arc` b ends around the Great Circle
            let b_gc_d = if Radians(0.0) < gc_d {
                b_d
            } else {
                core::f64::consts::TAU - b_d
            };
            if b_gc_d < abs_d.0 {
                // The end of the second `Arc` is beyond the end of first `Arc`
                (Radians(b_gc_d) + arc_length_0, arc_length_1)
            } else {
                // The start of the second `Arc` is before the start of first `Arc`
                (-abs_d, Radians(0.0))
            }
        }
    } else {
        // The distance to the start of arc2 from the end of arc1
        let b1a2 = if Radians(0.0) < gc_d {
            gc_d.0 - arc_length_0.0
        } else {
            core::f64::consts::TAU + gc_d.0 - arc_length_0.0
        };
        // The distance to the start of arc1 from the end of arc2
        let b2a1 = if Radians(0.0) < gc_d {
            core::f64::consts::TAU - gc_d.0 - arc_length_1.0
        } else {
            -gc_d.0 - arc_length_1.0
        };
        if b2a1 < b1a2 {
            // The start of the first arc is within the second arc
            (Radians(0.0), Radians(b2a1 + arc_length_1.0))
        } else {
            // The start of the second arc relative to the start of first arc.
            (Radians(b1a2 + arc_length_0.0), Radians(0.0))
        }
    }
}

/// Determine whether the antipodal point is closer to the centroid of the
/// `Arc`s.
///
/// * `point` a great-circle intersection point.
/// * `centroid` the centroid (geometric mean) of the `Arc`s mid points.
///
/// returns true if the antipodal intersection is closer to the `centroid`
/// of the `Arc`s otherwise returns false.
#[must_use]
pub fn use_antipodal_point(point: &Vector3d, centroid: &Vector3d) -> bool {
    sq_distance(centroid, &(-*point)) < sq_distance(centroid, point)
}

/// Calculate the great-circle distances along a pair of arcs to their
/// closest intersection point or their coincident arc distances if the
/// `Arc`s are on coincident Great Circles.
///
/// * `a_0`, `a_1` the `Arc` start points.
/// * `pole_0`, `pole_0` the `Arc` poles.
/// * `arc_length_0`, `arc_length_1` the `Arc` lengths.
/// * `centroid` the centroid (geometric mean) of the `Arc`s mid points.
///
/// returns the distances along the first arc and second arc to the intersection
/// point or to their coincident arc distances if the arcs do not intersect.
#[must_use]
pub fn calculate_intersection_point_distances(
    a_0: &Vector3d,
    pole_0: &Vector3d,
    arc_length_0: Radians,
    a_1: &Vector3d,
    pole_1: &Vector3d,
    arc_length_1: Radians,
    centroid: &Vector3d,
) -> (Radians, Radians) {
    // Calculate the square of the Euclidean distance between the start points.
    let sq_d = sq_distance(a_0, a_1);
    if sq_d < MIN_SQ_DISTANCE {
        (Radians(0.0), Radians(0.0))
    } else {
        calculate_intersection(pole_0, pole_1, MIN_SQ_NORM).map_or_else(
            || {
                calculate_coincident_arc_distances(
                    calculate_great_circle_atd(a_0, pole_0, a_1),
                    pole_0.dot(pole_1) < 0.0,
                    arc_length_0,
                    arc_length_1,
                )
            },
            |c| {
                // Find the closest intersection point
                let c = if use_antipodal_point(&c, centroid) {
                    -c
                } else {
                    c
                };
                calculate_intersection_distances(a_0, pole_0, a_1, pole_1, &c)
            },
        )
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
    // calculate the intersection point between the great circles
    let centroid = mid_point_0 + mid_point_1;
    let point = pole_0.cross(pole_1);
    normalise(&point, MIN_SQ_NORM).map_or_else(
        || {
            // the great circles are coincident

            let c = normalise_centroid(&centroid, mid_point_0, pole_0);

            // determine whether the great cicles oppose each other
            let angle = if pole_0.dot(pole_1).is_sign_negative() {
                Angle::default().opposite()
            } else {
                Angle::default()
            };

            (c, angle)
        },
        |p| {
            // find the closest intersection point to both the arcs
            let x = if use_antipodal_point(&p, &centroid) {
                -p
            } else {
                p
            };

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
    use crate::{LatLong, great_circle, vector};
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
    fn test_calculate_intersection_distances() {
        let start1 = LatLong::new(Degrees(-1.0), Degrees(-1.0));
        let a_0 = Vector3d::from(&start1);
        let azimuth1 = Angle::from(Degrees(45.0));
        let pole_0 = vector::calculate_pole(
            Angle::from(start1.lat()),
            Angle::from(start1.lon()),
            azimuth1,
        );

        let start2 = LatLong::new(Degrees(1.0), Degrees(-1.0));
        let a_1 = Vector3d::from(&start2);
        let azimuth2 = Angle::from(Degrees(135.0));
        let pole_1 = vector::calculate_pole(
            Angle::from(start2.lat()),
            Angle::from(start2.lon()),
            azimuth2,
        );

        let c = calculate_intersection(&pole_0, &pole_1, MIN_SQ_NORM).unwrap();
        let (c1, c2) = calculate_intersection_distances(&a_0, &pole_0, &a_1, &pole_1, &c);
        assert!(is_within_tolerance(-3.1169124762478333, c1.0, f64::EPSILON));
        assert!(is_within_tolerance(-3.1169124762478333, c2.0, f64::EPSILON));

        // Calculate the centre of the arc start points
        let centre_point = vector::normalise(&(a_0 + a_1), MIN_SQ_NORM).unwrap();
        assert!(sq_distance(&c, &centre_point) > 2.0);

        // opposite intersection point
        let d = -c;
        assert!(sq_distance(&d, &centre_point) <= 2.0);

        let (d1, d2) = calculate_intersection_distances(&a_0, &pole_0, &a_1, &pole_1, &d);
        assert!(is_within_tolerance(
            0.024680177341956263,
            d1.0,
            f64::EPSILON
        ));
        assert!(is_within_tolerance(
            0.024680177341956263,
            d2.0,
            f64::EPSILON
        ));

        // Same start points and intersection point
        let (e1, e2) = calculate_intersection_distances(&a_0, &pole_0, &a_0, &pole_1, &a_0);
        assert_eq!(0.0, e1.0);
        assert_eq!(0.0, e2.0);
    }

    #[test]
    fn test_is_alongside() {
        assert!(!is_alongside(
            Radians(-5.0 * f64::EPSILON),
            Radians(3.0),
            Radians(4.0 * f64::EPSILON)
        ));
        assert!(is_alongside(
            Radians(-4.0 * f64::EPSILON),
            Radians(3.0),
            Radians(4.0 * f64::EPSILON)
        ));
        assert!(is_alongside(
            Radians(3.0 + 4.0 * f64::EPSILON),
            Radians(3.0),
            Radians(4.0 * f64::EPSILON)
        ));
        assert!(!is_alongside(
            Radians(3.0 + 6.0 + f64::EPSILON),
            Radians(3.0),
            Radians(4.0 * f64::EPSILON)
        ));
    }

    #[test]
    fn test_is_within() {
        assert!(!is_within(-2.0 * f64::EPSILON, 2.0));
        assert!(is_within(-f64::EPSILON, 2.0));
        assert!(is_within(2.0 * (1.0 + f64::EPSILON), 2.0));
        assert!(!is_within(2.0 * (1.0 + 3.0 * f64::EPSILON), 2.0));
    }

    #[test]
    fn test_calculate_coincident_arc_distances() {
        let zero = Radians(0.0);
        let arc_length_0 = Radians(0.25);
        let arc_length_1 = Radians(0.75);

        let result0 =
            calculate_coincident_arc_distances(arc_length_1, true, arc_length_1, arc_length_0);
        assert_eq!(arc_length_1, result0.0);
        assert_eq!(zero, result0.1);

        let result1 =
            calculate_coincident_arc_distances(arc_length_1, true, arc_length_0, arc_length_1);
        assert_eq!(zero, result1.0);
        assert_eq!(arc_length_1, result1.1);

        let result2 =
            calculate_coincident_arc_distances(Radians(1.0), true, arc_length_0, arc_length_1);
        assert_eq!(arc_length_0, result2.0);
        assert_eq!(arc_length_1, result2.1);

        let result3 =
            calculate_coincident_arc_distances(Radians(1.5), true, arc_length_0, arc_length_1);
        assert_eq!(arc_length_1, result3.0);
        assert_eq!(arc_length_1, result3.1);

        let result4 =
            calculate_coincident_arc_distances(Radians(-1.5), true, arc_length_0, arc_length_1);
        assert_eq!(Radians(-1.5), result4.0);
        assert_eq!(zero, result4.1);

        let result5 =
            calculate_coincident_arc_distances(Radians(-1.0), false, arc_length_0, arc_length_1);
        assert_eq!(zero, result5.0);
        assert_eq!(Radians(1.0), result5.1);

        let result6 =
            calculate_coincident_arc_distances(Radians(1.0), false, arc_length_0, arc_length_1);
        assert_eq!(Radians(1.0), result6.0);
        assert_eq!(zero, result6.1);

        let result7 =
            calculate_coincident_arc_distances(-arc_length_1, false, arc_length_0, arc_length_1);
        assert_eq!(zero, result7.0);
        assert_eq!(arc_length_1, result7.1);

        let result8 =
            calculate_coincident_arc_distances(arc_length_0, false, arc_length_0, arc_length_1);
        assert_eq!(arc_length_0, result8.0);
        assert_eq!(zero, result8.1);
    }

    #[test]
    fn test_calculate_coincident_intersection_point_distances() {
        let start1 = LatLong::new(Degrees(0.0), Degrees(0.0));
        let a_0 = Vector3d::from(&start1);
        let end1 = LatLong::new(Degrees(0.0), Degrees(-4.0));
        let b1 = Vector3d::from(&end1);
        let pole_0 = &a_0.cross(&b1).normalize();
        let arc_length_0 = great_circle::e2gc_distance(vector::distance(&a_0, &b1));
        let mid_point1 = (a_0 + b1).normalize();

        let start2 = LatLong::new(Degrees(0.0), Degrees(0.25));
        let a_1 = Vector3d::from(&start2);
        let end2 = LatLong::new(Degrees(0.0), Degrees(4.0));
        let b2 = Vector3d::from(&end2);
        let pole_1 = &a_1.cross(&b2).normalize();
        let arc_length_1 = great_circle::e2gc_distance(vector::distance(&a_1, &b2));
        let mid_point2 = (a_1 + b2).normalize();

        let centriod = 0.5 * (mid_point1 + mid_point2);

        let distances = calculate_intersection_point_distances(
            &a_0,
            &pole_0,
            arc_length_0,
            &a_1,
            &pole_1,
            arc_length_1,
            &centriod,
        );
        assert!(is_within_tolerance(
            -0.25_f64.to_radians(),
            distances.0.0,
            f64::EPSILON
        ));
        assert_eq!(0.0, distances.1.0);

        let pole1r = -pole_0;
        let pole2r = -pole_1;

        let distances = calculate_intersection_point_distances(
            &b1,
            &pole1r,
            arc_length_0,
            &b2,
            &pole2r,
            arc_length_1,
            &centriod,
        );

        assert!(is_within_tolerance(
            4.25_f64.to_radians(),
            distances.0.0,
            f64::EPSILON
        ));
        assert_eq!(arc_length_1.0, distances.1.0);
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
