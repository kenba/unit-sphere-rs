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

//! # unit-sphere
//!
//! [![crates.io](https://img.shields.io/crates/v/unit-sphere.svg)](https://crates.io/crates/unit-sphere)
//! [![docs.io](https://docs.rs/unit-sphere/badge.svg)](https://docs.rs/unit-sphere/)
//! [![License](https://img.shields.io/badge/License-MIT-blue)](https://opensource.org/license/mit/)
//! [![Rust](https://github.com/kenba/unit-sphere-rs/actions/workflows/rust.yml/badge.svg)](https://github.com/kenba/unit-sphere-rs/actions)
//! [![codecov](https://codecov.io/gh/kenba/unit-sphere-rs/graph/badge.svg?token=G1H1XINERW)](https://codecov.io/gh/kenba/unit-sphere-rs)
//!
//! A library for performing geometric calculations on the surface of a sphere.
//!
//! The library uses a combination of spherical trigonometry and vector geometry
//! to perform [great-circle navigation](https://en.wikipedia.org/wiki/Great-circle_navigation)
//! on the surface of a unit sphere, see *Figure 1*.
//!
//! ![great circle path](https://upload.wikimedia.org/wikipedia/commons/thumb/c/cb/Illustration_of_great-circle_distance.svg/220px-Illustration_of_great-circle_distance.svg.png)  
//! *Figure 1 A Great Circle Path*
//!
//! A [great circle](https://en.wikipedia.org/wiki/Great_circle) is the
//! shortest path between positions on the surface of a sphere.  
//! It is the spherical equivalent of a straight line in planar geometry.
//!
//! ## Spherical trigonometry
//!
//! A great circle path between positions may be found using
//! [spherical trigonometry](https://en.wikipedia.org/wiki/Spherical_trigonometry).
//!
//! The [course](https://en.wikipedia.org/wiki/Great-circle_navigation#Course)
//! (initial azimuth) of a great circle can be calculated from the
//! latitudes and longitudes of the start and end points.  
//! While great circle distance can also be calculated from the latitudes and
//! longitudes of the start and end points using the
//! [haversine formula](https://en.wikipedia.org/wiki/Haversine_formula).  
//! The resulting distance in `Radians` can be converted to the required units by
//! multiplying the distance by the Earth radius measured in the required units.
//!
//! ## Vector geometry
//!
//! Points on the surface of a sphere and great circle poles may be represented
//! by 3D [vectors](https://www.movable-type.co.uk/scripts/latlong-vectors.html).  
//! Many calculations are simpler using vectors than spherical trigonometry.
//!
//! For example, the across track distance of a point from a great circle can
//! be calculated from the [dot product](https://en.wikipedia.org/wiki/Dot_product)
//! of the point and the great circle pole vectors.  
//! While intersection points of great circles can simply be calculated from
//! the [cross product](https://en.wikipedia.org/wiki/Cross_product) of their
//! pole vectors.
//!
//! ## Design
//!
//! The `great_circle` module performs spherical trigonometric calculations
//! and the `vector` module performs vector geometry calculations.
//!
//! The library is declared [no_std](https://docs.rust-embedded.org/book/intro/no-std.html)
//! so it can be used in embedded applications.

#![cfg_attr(not(test), no_std)]

extern crate angle_sc;
extern crate nalgebra as na;

pub mod great_circle;
pub mod vector;

use angle_sc::trig;
pub use angle_sc::{Angle, Degrees, Radians, Validate};

/// Test whether a latitude in degrees is a valid latitude.  
/// I.e. whether it lies in the range: -90.0 <= degrees <= 90.0
#[must_use]
pub fn is_valid_latitude(degrees: f64) -> bool {
    (-90.0..=90.0).contains(&degrees)
}

/// Test whether a longitude in degrees is a valid longitude.  
/// I.e. whether it lies in the range: -180.0 <= degrees <= 180.0
#[must_use]
pub fn is_valid_longitude(degrees: f64) -> bool {
    (-180.0..=180.0).contains(&degrees)
}

/// A position as a latitude and longitude pair of `Degrees`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LatLong {
    lat: Degrees,
    lon: Degrees,
}

impl Validate for LatLong {
    /// Test whether a `LatLong` is valid.  
    /// I.e. whether the latitude lies in the range: -90.0 <= lat <= 90.0
    /// and the longitude lies in the range: -90.0 <= lon <= 90.0
    #[must_use]
    fn is_valid(&self) -> bool {
        is_valid_latitude(self.lat.0) && is_valid_longitude(self.lon.0)
    }
}

impl LatLong {
    #[must_use]
    pub const fn new(lat: Degrees, lon: Degrees) -> Self {
        Self { lat, lon }
    }

    #[must_use]
    pub const fn lat(&self) -> Degrees {
        self.lat
    }

    #[must_use]
    pub const fn lon(&self) -> Degrees {
        self.lon
    }
}

impl TryFrom<(f64, f64)> for LatLong {
    type Error = &'static str;

    /// Attempt to convert a pair of f64 values in latitude, longitude order.
    ///
    /// return a valid `LatLong`.
    fn try_from(lat_long: (f64, f64)) -> Result<Self, Self::Error> {
        if !is_valid_latitude(lat_long.0) {
            Err("invalid latitude")
        } else if !is_valid_longitude(lat_long.1) {
            Err("invalid longitude")
        } else {
            Ok(Self::new(Degrees(lat_long.0), Degrees(lat_long.1)))
        }
    }
}

/// Calculate the azimuth and distance along the great circle of point b from
/// point a.
/// * `a`, `b` - the start and end positions
///
/// returns the Great Circle azimuth relative to North and distance of point b
/// from point a.
#[must_use]
pub fn calculate_azimuth_and_distance(a: &LatLong, b: &LatLong) -> (Angle, Radians) {
    let a_lat = Angle::from(a.lat);
    let b_lat = Angle::from(b.lat);
    let delta_long = Angle::from(b.lon - a.lon);
    (
        great_circle::calculate_gc_azimuth(a_lat, b_lat, delta_long),
        great_circle::calculate_gc_distance(a_lat, b_lat, delta_long),
    )
}

/// A `Vector3d` is a nalgebra Vector3.
#[allow(clippy::module_name_repetitions)]
pub type Vector3d = na::Vector3<f64>;

impl From<&LatLong> for Vector3d {
    /// Convert a `LatLong` to a point on the unit sphere
    /// @pre |lat| <= 90.0 degrees.
    /// * `lat` - the latitude.
    /// * `lon` - the longitude.
    ///
    /// returns a `Vector3d` of the point on the unit sphere.
    #[must_use]
    fn from(a: &LatLong) -> Self {
        vector::to_point(Angle::from(a.lat), Angle::from(a.lon))
    }
}

/// Calculate the latitude of a Point.
#[must_use]
pub fn latitude(a: &Vector3d) -> Angle {
    let sin_a = trig::UnitNegRange(a.z);
    Angle::new(sin_a, trig::swap_sin_cos(sin_a))
}

/// Calculate the longitude of a Point.
#[must_use]
pub fn longitude(a: &Vector3d) -> Angle {
    Angle::from_y_x(a.y, a.x)
}

impl From<&Vector3d> for LatLong {
    /// Convert a Point to a `LatLong`  
    #[must_use]
    fn from(value: &Vector3d) -> Self {
        Self::new(
            Degrees::from(latitude(value)),
            Degrees::from(longitude(value)),
        )
    }
}

/// An arc of a Great Circle on a unit sphere.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Arc {
    /// The start point of the arc.
    a: Vector3d,
    /// The right hand pole of the Great Circle of the arc.
    pole: Vector3d,
    /// The length of the arc.
    length: Radians,
    /// The half width of the arc.
    half_width: Radians,
}

impl Validate for Arc {
    /// Test whether an Arc is valid.
    /// I.e. both a and pole are on the unit sphere and are orthogonal and
    /// both length and `half_width` are >= 0.0.
    fn is_valid(&self) -> bool {
        vector::is_unit(&self.a)
            && vector::is_unit(&self.pole)
            && vector::are_orthogonal(&self.a, &self.pole)
            && (0.0 <= self.length.0)
            && (0.0 <= self.half_width.0)
    }
}

impl Arc {
    /// Construct an Arc
    /// * `a` - the start point of the arc.
    /// * `pole` - the right hand pole of the Great Circle of the arc.
    /// * `length` - the length of the arc.
    /// * `half_width` - the half width of the arc.
    #[must_use]
    pub const fn new(a: Vector3d, pole: Vector3d, length: Radians, half_width: Radians) -> Self {
        Self {
            a,
            pole,
            length,
            half_width,
        }
    }

    /// Construct an Arc
    /// * `a` - the start position
    /// * `azimuth` - the azimuth at a.
    /// * `length` - the length of the arc.
    #[must_use]
    pub fn from_lat_lon_azi_length(a: &LatLong, azimuth: Angle, length: Radians) -> Self {
        Self::new(
            Vector3d::from(a),
            vector::calculate_pole(Angle::from(a.lat()), Angle::from(a.lon()), azimuth),
            length,
            Radians(0.0),
        )
    }

    /// Construct an Arc from the start and end positions.  
    /// Note: if the points are the same or antipodal, the pole will be invalid.
    /// * `a`, `b` - the start and end positions
    #[must_use]
    pub fn between_positions(a: &LatLong, b: &LatLong) -> Self {
        let (azimuth, length) = calculate_azimuth_and_distance(a, b);
        let a_lat = Angle::from(a.lat());
        // if a is at the North or South pole
        if a_lat.cos().0 < great_circle::MIN_VALUE {
            // use b's longitude
            Self::from_lat_lon_azi_length(&LatLong::new(a.lat(), b.lon()), azimuth, length)
        } else {
            Self::from_lat_lon_azi_length(a, azimuth, length)
        }
    }

    /// Set the `half_width` of an Arc
    /// * `half_width` - the half width of the arc.
    #[must_use]
    pub fn set_half_width(&mut self, half_width: Radians) -> &mut Self {
        self.half_width = half_width;
        self
    }

    /// The start point of the arc.
    #[must_use]
    pub const fn a(&self) -> Vector3d {
        self.a
    }

    /// The right hand pole of the Great Circle at the start point of the arc.
    #[must_use]
    pub const fn pole(&self) -> Vector3d {
        self.pole
    }

    /// The length of the arc.
    #[must_use]
    pub const fn length(&self) -> Radians {
        self.length
    }

    /// The half width of the arc.
    #[must_use]
    pub const fn half_width(&self) -> Radians {
        self.half_width
    }

    /// The azimuth at the start point.      
    #[must_use]
    pub fn azimuth(&self) -> Angle {
        vector::calculate_azimuth(&self.a, &self.pole)
    }

    /// The direction vector of the arc at the start point.
    #[must_use]
    pub fn direction(&self) -> Vector3d {
        vector::direction(&self.a, &self.pole)
    }

    #[must_use]
    pub fn position(&self, distance: Radians) -> Vector3d {
        vector::position(&self.a, &self.direction(), Angle::from(distance))
    }

    /// The end point of the arc.
    #[must_use]
    pub fn b(&self) -> Vector3d {
        self.position(self.length)
    }

    /// The mid point of the arc.
    #[must_use]
    pub fn mid_point(&self) -> Vector3d {
        self.position(Radians(0.5 * self.length.0))
    }

    /// The position of a perpendicular point at distance from the arc.
    /// * `point` a point on the arc's great circle.
    /// * `distance` the perpendicular distance from the arc's great circle.
    ///
    /// returns the point at perpendicular distance from point.
    #[must_use]
    pub fn perp_position(&self, point: &Vector3d, distance: Radians) -> Vector3d {
        vector::position(point, &self.pole, Angle::from(distance))
    }

    /// The position of a point at angle from the arc start, at arc length.
    /// * `angle` the angle from the arc start.
    ///
    /// returns the point at angle from the arc start, at arc length.
    #[must_use]
    pub fn angle_position(&self, angle: Angle) -> Vector3d {
        vector::rotate_position(&self.a, &self.pole, angle, Angle::from(self.length))
    }

    /// The Arc at the end of an Arc, just the point if `half_width` is zero.
    /// @param `at_b` if true the arc at b, else the arc at a.
    ///
    /// @return the end arc at a or b.
    #[must_use]
    pub fn end_arc(&self, at_b: bool) -> Self {
        let p = if at_b { self.b() } else { self.a };
        let pole = vector::direction(&p, &self.pole);
        if self.half_width.0 < great_circle::MIN_VALUE {
            Self::new(p, pole, Radians(0.0), Radians(0.0))
        } else {
            let a = self.perp_position(&p, self.half_width);
            Self::new(a, pole, self.half_width + self.half_width, Radians(0.0))
        }
    }

    /// Calculate Great Circle along and across track distances of point from
    /// the Arc.
    /// * `point` - the point.
    ///
    /// returns the along and across track distances of the point in Radians.
    #[must_use]
    pub fn calculate_atd_and_xtd(&self, point: &Vector3d) -> (Radians, Radians) {
        vector::calculate_atd_and_xtd(&self.a, &self.pole(), point)
    }
}

impl TryFrom<(&LatLong, &LatLong)> for Arc {
    type Error = &'static str;

    /// Construct an Arc from a pair of positions.  
    /// * `params` - the start and end positions
    fn try_from(params: (&LatLong, &LatLong)) -> Result<Self, Self::Error> {
        // Convert positions to vectors
        let a = Vector3d::from(params.0);
        let b = Vector3d::from(params.1);
        // Calculate the great circle pole
        vector::normalize(&a.cross(&b)).map_or_else(
            || {
                if vector::sq_distance(&a, &b) < 1.0 {
                    Err("Positions are too close")
                } else {
                    Err("Positions are antipodal")
                }
            },
            |pole| {
                Ok(Self::new(
                    a,
                    pole,
                    great_circle::e2gc_distance(vector::distance(&a, &b)),
                    Radians(0.0),
                ))
            },
        )
    }
}

/// Calculate the distances along a pair of Arcs on the same (or reciprocal)
/// Great Circles to their closest intersection or reference points.
/// * `arc1`, `arc2` the arcs.
///
/// returns the distances along the first arc and second arc to the intersection
/// point or to their closest (reference) points if the arcs do not intersect.
#[must_use]
pub fn calculate_intersection_distances(arc1: &Arc, arc2: &Arc) -> (Radians, Radians) {
    vector::intersection::calculate_intersection_point_distances(
        &arc1.a,
        &arc1.pole,
        arc1.length(),
        &arc2.a,
        &arc2.pole,
        arc2.length(),
    )
}

/// Calculate whether a pair of Arcs intersect and (if so) where.
/// * `arc1`, `arc2` the arcs.
///
/// returns the distance along the first arc to the second arc or None if they
/// don't intersect.
#[must_use]
pub fn calculate_intersection_point_distance(arc1: &Arc, arc2: &Arc) -> Option<Radians> {
    let (distance1, distance2) = calculate_intersection_distances(arc1, arc2);

    // Determine whether both distances are within both arcs
    if vector::intersection::is_within(distance1.0, arc1.length().0)
        && vector::intersection::is_within(distance2.0, arc2.length().0)
    {
        Some(distance1)
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use angle_sc::{is_within_tolerance, Degrees};

    #[test]
    fn test_is_valid_latitude() {
        // value < -90
        assert!(!is_valid_latitude(-90.0001));
        // value = -90
        assert!(is_valid_latitude(-90.0));
        // value = 90
        assert!(is_valid_latitude(90.0));
        // value > 90
        assert!(!is_valid_latitude(90.0001));
    }

    #[test]
    fn test_is_valid_longitude() {
        // value < -180
        assert!(!is_valid_longitude(-180.0001));
        // value = -180
        assert!(is_valid_longitude(-180.0));
        // value = 180
        assert!(is_valid_longitude(180.0));
        // value > 180
        assert!(!is_valid_longitude(180.0001));
    }

    #[test]
    fn test_latlong_traits() {
        let a = LatLong::try_from((0.0, 90.0)).unwrap();

        assert!(a.is_valid());

        let a_clone = a.clone();
        assert!(a_clone == a);

        assert_eq!(Degrees(0.0), a.lat());
        assert_eq!(Degrees(90.0), a.lon());

        print!("LatLong: {:?}", a);

        let invalid_lat = LatLong::try_from((91.0, 0.0));
        assert_eq!(Err("invalid latitude"), invalid_lat);

        let invalid_lon = LatLong::try_from((0.0, 181.0));
        assert_eq!(Err("invalid longitude"), invalid_lon);
    }

    #[test]
    fn test_vector3d_traits() {
        let a = LatLong::try_from((0.0, 90.0)).unwrap();
        let point = Vector3d::from(&a);

        assert_eq!(0.0, point.x);
        assert_eq!(1.0, point.y);
        assert_eq!(0.0, point.z);

        assert_eq!(Degrees(0.0), Degrees::from(latitude(&point)));
        assert_eq!(Degrees(90.0), Degrees::from(longitude(&point)));

        let result = LatLong::from(&point);
        assert_eq!(a, result);
    }

    #[test]
    fn test_great_circle_90n_0n_0e() {
        let a = LatLong::new(Degrees(90.0), Degrees(0.0));
        let b = LatLong::new(Degrees(0.0), Degrees(0.0));
        let (azimuth, dist) = calculate_azimuth_and_distance(&a, &b);

        assert!(is_within_tolerance(
            core::f64::consts::FRAC_PI_2,
            dist.0,
            f64::EPSILON
        ));
        assert_eq!(180.0, Degrees::from(azimuth).0);
    }

    #[test]
    fn test_great_circle_90s_0n_50e() {
        let a = LatLong::new(Degrees(-90.0), Degrees(0.0));
        let b = LatLong::new(Degrees(0.0), Degrees(50.0));
        let (azimuth, dist) = calculate_azimuth_and_distance(&a, &b);

        assert!(is_within_tolerance(
            core::f64::consts::FRAC_PI_2,
            dist.0,
            f64::EPSILON
        ));
        assert_eq!(0.0, Degrees::from(azimuth).0);
    }

    #[test]
    fn test_great_circle_0n_60e_0n_60w() {
        let a = LatLong::new(Degrees(0.0), Degrees(60.0));
        let b = LatLong::new(Degrees(0.0), Degrees(-60.0));
        let (azimuth, dist) = calculate_azimuth_and_distance(&a, &b);

        assert!(is_within_tolerance(
            2.0 * core::f64::consts::FRAC_PI_3,
            dist.0,
            2.0 * f64::EPSILON
        ));
        assert_eq!(-90.0, Degrees::from(azimuth).0);
    }

    #[test]
    fn test_arc() {
        // Greenwich equator
        let g_eq = LatLong::new(Degrees(0.0), Degrees(0.0));

        // 90 degrees East on the equator
        let e_eq = LatLong::new(Degrees(0.0), Degrees(90.0));

        let mut arc = Arc::between_positions(&g_eq, &e_eq);
        let arc = arc.set_half_width(Radians(0.01));
        assert!(arc.is_valid());
        assert_eq!(Radians(0.01), arc.half_width());

        assert_eq!(Vector3d::from(&g_eq), arc.a());
        assert_eq!(Vector3d::new(0.0, 0.0, 1.0), arc.pole());
        assert!(is_within_tolerance(
            core::f64::consts::FRAC_PI_2,
            arc.length().0,
            f64::EPSILON
        ));
        assert_eq!(Angle::from(Degrees(90.0)), arc.azimuth());
        let b = Vector3d::from(&e_eq);
        assert!(is_within_tolerance(
            0.0,
            vector::distance(&b, &arc.b()),
            f64::EPSILON
        ));

        let mid_point = arc.mid_point();
        assert_eq!(0.0, mid_point.z);
        assert!(is_within_tolerance(
            45.0,
            Degrees::from(longitude(&mid_point)).0,
            32.0 * f64::EPSILON
        ));

        let start_arc = arc.end_arc(false);
        assert_eq!(0.02, start_arc.length().0);

        let start_arc_a = start_arc.a();
        assert_eq!(start_arc_a, arc.perp_position(&arc.a(), Radians(0.01)));

        let angle_90 = Angle::from(Degrees(90.0));
        let pole_0 = Vector3d::new(0.0, 0.0, 1.0);
        assert!(vector::distance(&pole_0, &arc.angle_position(angle_90)) <= f64::EPSILON);

        let end_arc = arc.end_arc(true);
        assert_eq!(0.02, end_arc.length().0);

        let end_arc_a = end_arc.a();
        assert_eq!(end_arc_a, arc.perp_position(&arc.b(), Radians(0.01)));
    }

    #[test]
    fn test_north_and_south_poles() {
        let north_pole = LatLong::new(Degrees(90.0), Degrees(0.0));
        let south_pole = LatLong::new(Degrees(-90.0), Degrees(0.0));

        let (azimuth, distance) = calculate_azimuth_and_distance(&south_pole, &north_pole);
        assert_eq!(0.0, Degrees::from(azimuth).0);
        assert_eq!(core::f64::consts::PI, distance.0);

        let (azimuth, distance) = calculate_azimuth_and_distance(&north_pole, &south_pole);
        assert_eq!(180.0, Degrees::from(azimuth).0);
        assert_eq!(core::f64::consts::PI, distance.0);

        // 90 degrees East on the equator
        let e_eq = LatLong::new(Degrees(0.0), Degrees(50.0));

        let arc = Arc::between_positions(&north_pole, &e_eq);
        assert!(is_within_tolerance(
            e_eq.lat().0,
            LatLong::from(&arc.b()).lat().abs().0,
            1e-13
        ));
        assert!(is_within_tolerance(
            e_eq.lon().0,
            LatLong::from(&arc.b()).lon().0,
            50.0 * f64::EPSILON
        ));

        let arc = Arc::between_positions(&south_pole, &e_eq);
        assert!(is_within_tolerance(
            e_eq.lat().0,
            LatLong::from(&arc.b()).lat().abs().0,
            1e-13
        ));
        assert!(is_within_tolerance(
            e_eq.lon().0,
            LatLong::from(&arc.b()).lon().0,
            50.0 * f64::EPSILON
        ));

        let w_eq = LatLong::new(Degrees(0.0), Degrees(-140.0));

        let arc = Arc::between_positions(&north_pole, &w_eq);
        assert!(is_within_tolerance(
            w_eq.lat().0,
            LatLong::from(&arc.b()).lat().abs().0,
            1e-13
        ));
        assert!(is_within_tolerance(
            w_eq.lon().0,
            LatLong::from(&arc.b()).lon().0,
            256.0 * f64::EPSILON
        ));

        let arc = Arc::between_positions(&south_pole, &w_eq);
        assert!(is_within_tolerance(
            w_eq.lat().0,
            LatLong::from(&arc.b()).lat().abs().0,
            1e-13
        ));
        assert!(is_within_tolerance(
            w_eq.lon().0,
            LatLong::from(&arc.b()).lon().0,
            256.0 * f64::EPSILON
        ));

        let invalid_arc = Arc::try_from((&north_pole, &north_pole));
        assert_eq!(Err("Positions are too close"), invalid_arc);

        let arc = Arc::between_positions(&north_pole, &north_pole);
        assert_eq!(north_pole, LatLong::from(&arc.b()));

        let invalid_arc = Arc::try_from((&north_pole, &south_pole));
        assert_eq!(Err("Positions are antipodal"), invalid_arc);

        let arc = Arc::between_positions(&north_pole, &south_pole);
        assert_eq!(south_pole, LatLong::from(&arc.b()));

        let arc = Arc::between_positions(&south_pole, &north_pole);
        assert_eq!(north_pole, LatLong::from(&arc.b()));

        let arc = Arc::between_positions(&south_pole, &south_pole);
        assert_eq!(south_pole, LatLong::from(&arc.b()));
    }

    #[test]
    fn test_arc_atd_and_xtd() {
        // Greenwich equator
        let g_eq = LatLong::new(Degrees(0.0), Degrees(0.0));

        // 90 degrees East on the equator
        let e_eq = LatLong::new(Degrees(0.0), Degrees(90.0));

        let arc = Arc::try_from((&g_eq, &e_eq)).unwrap();
        assert!(arc.is_valid());

        let start_arc = arc.end_arc(false);
        assert_eq!(0.0, start_arc.length().0);

        let start_arc_a = start_arc.a();
        assert_eq!(arc.a(), start_arc_a);

        let longitude = Degrees(1.0);

        // Test across track distance
        // Accuracy drops off outside of this range
        for lat in -83..84 {
            let latitude = Degrees(lat as f64);
            let latlong = LatLong::new(latitude, longitude);
            let point = Vector3d::from(&latlong);

            let expected = (lat as f64).to_radians();
            let (atd, xtd) = arc.calculate_atd_and_xtd(&point);
            assert!(is_within_tolerance(1_f64.to_radians(), atd.0, f64::EPSILON));
            assert!(is_within_tolerance(expected, xtd.0, 2.0 * f64::EPSILON));
        }
    }

    #[test]
    fn test_arc_intersection_point_length() {
        // Karney's example
        // Istanbul, Washington, Reyjavik and Accra
        let istanbul = LatLong::new(Degrees(42.0), Degrees(29.0));
        let washington = LatLong::new(Degrees(39.0), Degrees(-77.0));
        let reyjavik = LatLong::new(Degrees(64.0), Degrees(-22.0));
        let accra = LatLong::new(Degrees(6.0), Degrees(0.0));

        let arc1 = Arc::try_from((&istanbul, &washington)).unwrap();
        let arc2 = Arc::try_from((&reyjavik, &accra)).unwrap();

        let intersection_distance = calculate_intersection_point_distance(&arc1, &arc2).unwrap();
        assert!(is_within_tolerance(
            0.5406004765152588,
            intersection_distance.0,
            f64::EPSILON
        ));

        let intersection_distance_other_arc =
            calculate_intersection_point_distance(&arc2, &arc1).unwrap();
        assert!(is_within_tolerance(
            0.17553891720631054,
            intersection_distance_other_arc.0,
            f64::EPSILON
        ));

        let intersection_pos = arc1.position(intersection_distance);
        let lat_long = LatLong::from(&intersection_pos);
        // Geodesic intersection latitude is 54.7170296111
        assert!(is_within_tolerance(54.7, lat_long.lat().0, 0.4));
        // Geodesic intersection longitude is -14.56385575
        assert!(is_within_tolerance(-14.56, lat_long.lon().0, 0.02));
    }

    #[test]
    fn test_arc_intersection_same_gerat_circles() {
        let south_pole_1 = LatLong::new(Degrees(-88.0), Degrees(-180.0));
        let south_pole_2 = LatLong::new(Degrees(-87.0), Degrees(0.0));

        let arc1 = Arc::try_from((&south_pole_1, &south_pole_2)).unwrap();

        let intersection_lengths = calculate_intersection_distances(&arc1, &arc1);
        assert_eq!(Radians(0.0), intersection_lengths.0);
        assert_eq!(Radians(0.0), intersection_lengths.1);

        let gc_length = calculate_intersection_point_distance(&arc1, &arc1).unwrap();
        assert_eq!(Radians(0.0), gc_length);

        let south_pole_3 = LatLong::new(Degrees(-85.0), Degrees(0.0));
        let south_pole_4 = LatLong::new(Degrees(-86.0), Degrees(0.0));
        let arc2 = Arc::try_from((&south_pole_3, &south_pole_4)).unwrap();
        let intersection_length = calculate_intersection_point_distance(&arc1, &arc2);
        assert!(intersection_length.is_none());
    }
}
