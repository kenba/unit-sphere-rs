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

//! The `intersection` module contains functions for calculating great-circle
//! intersections using vectors.
//!
//! A pair of great circles intersect at two points unless they are coincident.  
//! For example, points `u` and `v` in *Figure1*.
//!
//! ![great circle path](https://upload.wikimedia.org/wikipedia/commons/thumb/c/cb/Illustration_of_great-circle_distance.svg/220px-Illustration_of_great-circle_distance.svg.png)  
//! *Figure 1 A pair of intersecting great circles*
//!
//! A great circle intersection point can simply be calculated by normalizing
//! the [cross product](https://en.wikipedia.org/wiki/Cross_product) of their
//! pole vectors.  
//! If the resulting vector is too small to normalize, then the great circles
//! are coincident, in which case they effectively *intersect* everywhere.
//!
//! If a pair of `Arc`s are on coincident great circles,
//! `calculate_coincident_arc_distances` calculates the distances between
//! `Arc` ends, zero if the `Arc`s overlap.
//! 
//! Otherwise `closest_intersection_point` determines which intersection point
//! is closer to the [centroid](https://en.wikipedia.org/wiki/Centroid)
//! of the `Arc`s midpoints, normalized to lie on the surface of the sphere.
//! `calculate_intersection_distances` then calculates great-circle distances
//! along the `Arc`s to the intersection point.

use super::{calculate_great_circle_atd, normalise, sq_distance, Vector3d, MIN_SQ_DISTANCE};
use angle_sc::{max, Radians};

/// Calculate an intersection point between the poles of two Great Circles.  
/// See: <http://www.movable-type.co.uk/scripts/latlong-vectors.html#intersection>  
/// * `pole1`, `pole2` the poles.
///
/// return an intersection point or None if the poles represent coincident Great Circles.
#[must_use]
pub fn calculate_intersection_point(pole1: &Vector3d, pole2: &Vector3d) -> Option<Vector3d> {
    normalise(&pole1.cross(pole2))
}

/// Calculate the great circle distances to an intersection point from the
/// start points of a pair of great circle arcs, on different great circles.
/// * `a1`, `a2` the start points of the great circle arcs
/// * `pole1`, `pole2` the poles of the great circle arcs
/// * `c` the intersection point
///
/// returns a pair of great circle distances along the arcs to the
/// intersection point in `Radians`.
#[must_use]
pub fn calculate_intersection_distances(
    a1: &Vector3d,
    pole1: &Vector3d,
    a2: &Vector3d,
    pole2: &Vector3d,
    c: &Vector3d,
) -> (Radians, Radians) {
    (
        calculate_great_circle_atd(a1, pole1, c),
        calculate_great_circle_atd(a2, pole2, c),
    )
}

/// Whether an intersection point is within an `Arc`.
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
/// * `gc_d` the great-circle distance between the arc start points.
/// * `reciprocal` whether the arcs are in reciprocal directions.
/// * `arc1_length`, `arc2_length` the `Arc` lengths in `Radians`.
///
/// returns the distances along the first `Arc` and second `Arc` to their closest
/// (reference) points in `Radians`.
#[must_use]
pub fn calculate_coincident_arc_distances(
    gc_d: Radians,
    reciprocal: bool,
    arc1_length: Radians,
    arc2_length: Radians,
) -> (Radians, Radians) {
    if reciprocal {
        // if the arcs intersect
        if is_within(gc_d.0, max(arc1_length, arc2_length).0) {
            if gc_d <= arc2_length {
                // The start of the first `Arc` is within the second `Arc`
                (Radians(0.0), gc_d)
            } else {
                // The start of the second `Arc` is within the first `Arc`
                (gc_d, Radians(0.0))
            }
        } else {
            let abs_d = gc_d.abs();

            // The distance between the `Arc` b ends
            let b_d = abs_d.0 - arc1_length.0 - arc2_length.0;
            // The distance between the `Arc` b ends around the Great Circle
            let b_gc_d = if Radians(0.0) < gc_d {
                b_d
            } else {
                core::f64::consts::TAU - b_d
            };
            if b_gc_d < abs_d.0 {
                // The end of the second `Arc` is beyond the end of first `Arc`
                (Radians(b_gc_d) + arc1_length, arc2_length)
            } else {
                // The start of the second `Arc` is before the start of first `Arc`
                (-abs_d, Radians(0.0))
            }
        }
    } else {
        // The distance to the start of arc2 from the end of arc1
        let b1a2 = if Radians(0.0) < gc_d {
            gc_d.0 - arc1_length.0
        } else {
            core::f64::consts::TAU + gc_d.0 - arc1_length.0
        };
        // The distance to the start of arc1 from the end of arc2
        let b2a1 = if Radians(0.0) < gc_d {
            core::f64::consts::TAU - gc_d.0 - arc2_length.0
        } else {
            -gc_d.0 - arc2_length.0
        };
        if b2a1 < b1a2 {
            // The start of the first arc is within the second arc
            (Radians(0.0), Radians(b2a1 + arc2_length.0))
        } else {
            // The start of the second arc relative to the start of first arc.
            (Radians(b1a2 + arc1_length.0), Radians(0.0))
        }
    }
}

/// Determine the closest intersection point to the centre of the `Arc`s.
/// * `intersection` a great circle intersection point.
/// * `centre` the centre of the `Arc`s mid points.
///
/// returns `intersection` if `centre` cannot be normalized or `intersection`
/// and `centre` are in the same hemisphere. Otherwise, it returns the
/// antipodal intersection point.
#[must_use]
pub fn closest_intersection_point(intersection: &Vector3d, centre: Option<Vector3d>) -> Vector3d {
    let use_antipodal_point = centre.is_some_and(|x| sq_distance(intersection, &x) > 2.0);
    if use_antipodal_point {
        -*intersection
    } else {
        *intersection
    }
}

/// Calculate the great-circle distances along a pair of arcs to their
/// closest intersection point or their coincident arc distances if the
/// `Arc`s are on coincident Great Circles.
/// * `a1`, `a2` the `Arc` start points.
/// * `pole1`, `pole1` the `Arc` poles.
/// * `length1`, `length2` the `Arc` lengths.
/// * `centre` the centre of the `Arc`s mid points.
///
/// returns the distances along the first arc and second arc to the intersection
/// point or to their coincident arc distances if the arcs do not intersect.
#[must_use]
pub fn calculate_intersection_point_distances(
    a1: &Vector3d,
    pole1: &Vector3d,
    length1: Radians,
    a2: &Vector3d,
    pole2: &Vector3d,
    length2: Radians,
    centre: Option<Vector3d>,
) -> (Radians, Radians) {
    // Calculate the square of the Euclidean distance between the start points.
    let sq_d = sq_distance(a1, a2);
    if sq_d < MIN_SQ_DISTANCE {
        (Radians(0.0), Radians(0.0))
    } else {
        calculate_intersection_point(pole1, pole2).map_or_else(
            || {
                calculate_coincident_arc_distances(
                    calculate_great_circle_atd(a1, pole1, a2),
                    pole1.dot(pole2) < 0.0,
                    length1,
                    length2,
                )
            },
            |c| {
                calculate_intersection_distances(
                    a1,
                    pole1,
                    a2,
                    pole2,
                    &closest_intersection_point(&c, centre),
                )
            },
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{vector, LatLong};
    use angle_sc::{is_within_tolerance, Angle, Degrees};

    #[test]
    fn test_calculate_intersection_point() {
        let lat_lon_south = LatLong::new(Degrees(-90.0), Degrees(0.0));
        let south_pole = Vector3d::from(&lat_lon_south);

        let lat_lon_north = LatLong::new(Degrees(90.0), Degrees(0.0));
        let north_pole = Vector3d::from(&lat_lon_north);

        let lat_lon_idl = LatLong::new(Degrees(0.0), Degrees(180.0));
        let idl = Vector3d::from(&lat_lon_idl);

        let equator_intersection = calculate_intersection_point(&south_pole, &north_pole);
        assert!(equator_intersection.is_none());

        let gc_intersection1 = calculate_intersection_point(&idl, &north_pole).unwrap();
        let gc_intersection2 = calculate_intersection_point(&idl, &south_pole).unwrap();

        assert_eq!(gc_intersection1, -gc_intersection2);
    }

    #[test]
    fn test_calculate_intersection_distances() {
        let start1 = LatLong::new(Degrees(-1.0), Degrees(-1.0));
        let a1 = Vector3d::from(&start1);
        let azimuth1 = Angle::from(Degrees(45.0));
        let pole1 = vector::calculate_pole(
            Angle::from(start1.lat()),
            Angle::from(start1.lon()),
            azimuth1,
        );

        let start2 = LatLong::new(Degrees(1.0), Degrees(-1.0));
        let a2 = Vector3d::from(&start2);
        let azimuth2 = Angle::from(Degrees(135.0));
        let pole2 = vector::calculate_pole(
            Angle::from(start2.lat()),
            Angle::from(start2.lon()),
            azimuth2,
        );

        let c = calculate_intersection_point(&pole1, &pole2).unwrap();
        let (c1, c2) = calculate_intersection_distances(&a1, &pole1, &a2, &pole2, &c);
        assert!(is_within_tolerance(-3.1169124762478333, c1.0, f64::EPSILON));
        assert!(is_within_tolerance(-3.1169124762478333, c2.0, f64::EPSILON));

        // Calculate the centre of the arc start points
        let centre_point = vector::normalise(&(a1 + a2)).unwrap();
        assert!(sq_distance(&c, &centre_point) > 2.0);

        // opposite intersection point
        let d = -c;
        assert!(sq_distance(&d, &centre_point) <= 2.0);

        let (d1, d2) = calculate_intersection_distances(&a1, &pole1, &a2, &pole2, &d);
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
        let (e1, e2) = calculate_intersection_distances(&a1, &pole1, &a1, &pole2, &a1);
        assert_eq!(0.0, e1.0);
        assert_eq!(0.0, e2.0);
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
        let length1 = Radians(0.25);
        let length2 = Radians(0.75);

        let result0 = calculate_coincident_arc_distances(length2, true, length2, length1);
        assert_eq!(length2, result0.0);
        assert_eq!(zero, result0.1);

        let result1 = calculate_coincident_arc_distances(length2, true, length1, length2);
        assert_eq!(zero, result1.0);
        assert_eq!(length2, result1.1);

        let result2 = calculate_coincident_arc_distances(Radians(1.0), true, length1, length2);
        assert_eq!(length1, result2.0);
        assert_eq!(length2, result2.1);

        let result3 = calculate_coincident_arc_distances(Radians(1.5), true, length1, length2);
        assert_eq!(length2, result3.0);
        assert_eq!(length2, result3.1);

        let result4 = calculate_coincident_arc_distances(Radians(-1.5), true, length1, length2);
        assert_eq!(Radians(-1.5), result4.0);
        assert_eq!(zero, result4.1);

        let result5 = calculate_coincident_arc_distances(Radians(-1.0), false, length1, length2);
        assert_eq!(zero, result5.0);
        assert_eq!(Radians(1.0), result5.1);

        let result6 = calculate_coincident_arc_distances(Radians(1.0), false, length1, length2);
        assert_eq!(Radians(1.0), result6.0);
        assert_eq!(zero, result6.1);

        let result7 = calculate_coincident_arc_distances(-length2, false, length1, length2);
        assert_eq!(zero, result7.0);
        assert_eq!(length2, result7.1);

        let result8 = calculate_coincident_arc_distances(length1, false, length1, length2);
        assert_eq!(length1, result8.0);
        assert_eq!(zero, result8.1);
    }
}
