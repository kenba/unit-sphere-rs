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

//! The `great_circle` module contains functions for calculating the course
//! and distance between points along great circles on a unit sphere.

use angle_sc::{trig, Angle, Radians};

/// The minimum value for angles and distances.
pub const MIN_VALUE: f64 = 2.0 * f64::EPSILON;

/// Calculate the Great Circle distance between two points from their
/// Latitude and Longitude differences.
///
/// See: [Haversine formula](https://en.wikipedia.org/wiki/Haversine_formula).
/// This function is less accurate than `calculate_gc_distance`.
/// * `a_lat` - start point Latitude.
/// * `b_lat` - finish point Latitude.
/// * `delta_long` - Longitude difference between start and finish points.
/// * `delta_lat` - Latitude difference between start and finish points.
///
/// returns the Great Circle distance between the points in Radians.
#[must_use]
pub fn calculate_haversine_distance(
    a_lat: Angle,
    b_lat: Angle,
    delta_long: Angle,
    delta_lat: Angle,
) -> Radians {
    let haversine_lat = trig::sq_sine_half(delta_lat.cos());
    let haversine_lon = trig::sq_sine_half(delta_long.cos());

    let a = (haversine_lat + a_lat.cos().0 * b_lat.cos().0 * haversine_lon).clamp(0.0, 1.0);

    if a < MIN_VALUE {
        Radians(0.0)
    } else {
        Radians(2.0 * libm::asin(libm::sqrt(a)))
    }
}

/// Convert a Euclidean distance to a Great Circle distance (in `Radians`).
/// e should satisfy: 0 <= e <= 2, if not it is clamped into range.
#[must_use]
pub fn e2gc_distance(e: f64) -> Radians {
    if e < MIN_VALUE {
        Radians(0.0)
    } else {
        Radians(2.0 * libm::asin(trig::UnitNegRange::clamp(0.5 * e).0))
    }
}

/// Convert a Great Circle distance (in radians) to a Euclidean distance.
#[must_use]
pub fn gc2e_distance(gc: Radians) -> f64 {
    2.0 * libm::sin(0.5 * gc.0)
}

/// Calculate the square of the Euclidean distance (i.e. using Pythagoras)
/// between two points from their Latitudes and their Longitude difference.
///
/// * `a_lat` - start point Latitude.
/// * `b_lat` - finish point Latitude.
/// * `delta_long` - Longitude difference between start and finish points.
///
/// returns the square of the Euclidean distance between the points.
#[must_use]
pub fn sq_euclidean_distance(a_lat: Angle, b_lat: Angle, delta_long: Angle) -> f64 {
    let delta_x = b_lat.cos().0 * delta_long.cos().0 - a_lat.cos().0;
    let delta_y = b_lat.cos().0 * delta_long.sin().0;
    let delta_z = b_lat.sin().0 - a_lat.sin().0;

    let result = delta_x * delta_x + delta_y * delta_y + delta_z * delta_z;
    result.clamp(0.0, 4.0)
}

/// Calculate the Great Circle distance (angle from centre) between two points
/// from their Latitudes and their Longitude difference.
///
/// This function is more accurate than `calculate_haversine_distance`.
/// * `a_lat` - start point Latitude.
/// * `b_lat` - finish point Latitude.
/// * `delta_long` - Longitude difference between start and finish points.
///
/// returns the Great Circle distance between the points in Radians.
#[must_use]
pub fn calculate_gc_distance(a_lat: Angle, b_lat: Angle, delta_long: Angle) -> Radians {
    e2gc_distance(libm::sqrt(sq_euclidean_distance(a_lat, b_lat, delta_long)))
}

/// Calculate the azimuth (bearing) along the great circle of point b from
/// point a from their Latitudes and their Longitude difference.
///
/// * `a_lat` - start point Latitude.
/// * `b_lat` - finish point Latitude.
/// * `delta_long` - Longitude difference between start and finish points.
///
/// returns the Great Circle azimuth relative to North of point b from point a
/// as an Angle.
#[must_use]
pub fn calculate_gc_azimuth(a_lat: Angle, b_lat: Angle, delta_long: Angle) -> Angle {
    // if start point is North or South pole
    if a_lat.cos().0 < MIN_VALUE {
        if a_lat.sin().0.is_sign_negative() {
            // South pole, azimuth is zero
            Angle::default()
        } else {
            // North pole, azimuth is 180 degrees
            Angle::new(trig::UnitNegRange(0.0), trig::UnitNegRange(-1.0))
        }
    } else {
        let sin_azimuth = b_lat.cos().0 * delta_long.sin().0;
        let temp = (a_lat.sin().0 * b_lat.cos().0 * delta_long.sin().0 * delta_long.sin().0)
            / (1.0 + delta_long.cos().abs().0);
        let cos_azimuth = if delta_long.cos().0 < 0.0 {
            b_lat.sin().0 * a_lat.cos().0 + a_lat.sin().0 * b_lat.cos().0 - temp
        } else {
            b_lat.sin().0 * a_lat.cos().0 - a_lat.sin().0 * b_lat.cos().0 + temp
        };

        Angle::from_y_x(sin_azimuth, cos_azimuth)
    }
}

/// Calculate the Great Circle distance (as an angle, sigma) between two points
/// from their Latitudes and their Longitude difference.
///
/// From Eurocae ED-323 Appendix E.
/// * `a_lat` - start point latitude.
/// * `b_lat` - finish point latitude.
/// * `delta_long` - longitude difference between start and finish points.
///
/// returns the Great Circle distance between the points (sigma) as an Angle.
pub fn calculate_sigma(a_lat: Angle, b_lat: Angle, delta_long: Angle) -> Angle {
    let a = b_lat.cos().0 * delta_long.sin().0;
    let b = a_lat.cos().0 * b_lat.sin().0 - a_lat.sin().0 * b_lat.cos().0 * delta_long.cos().0;
    let sin_sigma = trig::UnitNegRange(libm::hypot(a, b));
    let cos_sigma = trig::UnitNegRange(
        a_lat.sin().0 * b_lat.sin().0 + a_lat.cos().0 * b_lat.cos().0 * delta_long.cos().0,
    );
    Angle::new(sin_sigma, cos_sigma)
}

/// Calculate the latitude at great circle distance, sigma.
/// * `a_lat` - the latitude of the start point.
/// * `azi` - the azimuth at a_lat.
/// * `sigma` - the distance on the auxiliary sphere as an Angle.
///
/// return the latitude at sigma from a_lat.
#[must_use]
pub fn calculate_latitude(lat: Angle, azi: Angle, sigma: Angle) -> Angle {
    let sin_lat =
        trig::UnitNegRange(lat.sin().0 * sigma.cos().0 + lat.cos().0 * sigma.sin().0 * azi.cos().0);
    let cos_lat = trig::UnitNegRange(libm::hypot(
        lat.cos().0 * azi.sin().0,
        lat.sin().0 * sigma.sin().0 - lat.cos().0 * sigma.cos().0 * azi.cos().0,
    ));

    Angle::new(sin_lat, cos_lat)
}

/// Calculate the longitude difference at great circle distance, sigma.
/// * `a_lat` - the latitude of the start point.
/// * `azi` - the azimuth at a_lat.
/// * `sigma` - the distance on the auxiliary sphere as an Angle.
///
/// return the longitude difference from a_lat.
#[must_use]
pub fn calculate_delta_longitude(lat: Angle, azi: Angle, sigma: Angle) -> Angle {
    let sin_lon = trig::UnitNegRange(sigma.sin().0 * azi.sin().0);
    let cos_lon =
        trig::UnitNegRange(lat.cos().0 * sigma.cos().0 - lat.sin().0 * sigma.sin().0 * azi.cos().0);

    Angle::new(sin_lon, cos_lon)
}

/// Calculate the azimuth at latitude b_lat given the latitude a_lat and azimuth, `alpha1`.
/// * `a_lat`, `b_lat` - the latitudes of the points on the unit sphere.
/// * `azi` - the azimuth at a_lat.
///
/// returns the azimuth at b_lat.
#[must_use]
pub fn calculate_other_azimuth(a_lat: Angle, b_lat: Angle, azi: Angle) -> Angle {
    let clairaut = trig::UnitNegRange(azi.sin().0 * a_lat.cos().0);

    let sin_alpha2 = if b_lat.cos() == a_lat.cos() {
        azi.sin()
    } else {
        trig::UnitNegRange::clamp(clairaut.0 / b_lat.cos().0)
    };

    // Karney's method to calculate the cosine of the other azimuth
    let cos_alpha2 = if (b_lat.cos() != a_lat.cos()) || (b_lat.sin().abs().0 != -a_lat.sin().0) {
        let temp1 = azi.cos().0 * a_lat.cos().0;
        let temp2 = if a_lat.cos().0 < a_lat.abs().sin().0 {
            trig::sq_a_minus_sq_b(b_lat.cos(), a_lat.cos()).0
        } else {
            trig::sq_a_minus_sq_b(a_lat.sin(), b_lat.sin()).0
        };
        let temp3 = temp1 * temp1 + temp2;
        let temp4 = if temp3 > 0.0 {
            libm::sqrt(temp3) / b_lat.cos().0
        } else {
            0.0
        };
        trig::UnitNegRange::clamp(temp4)
    } else {
        azi.cos().abs()
    };

    Angle::new(sin_alpha2, cos_alpha2)
}

#[cfg(test)]
mod tests {
    use super::*;
    use angle_sc::{is_within_tolerance, Degrees};

    #[test]
    fn test_distance_conversion_functions() {
        assert_eq!(core::f64::consts::PI, e2gc_distance(2.0).0);
        assert!(is_within_tolerance(
            core::f64::consts::FRAC_PI_2,
            e2gc_distance(core::f64::consts::SQRT_2).0,
            f64::EPSILON
        ));
        assert!(is_within_tolerance(
            core::f64::consts::SQRT_2,
            gc2e_distance(Radians(core::f64::consts::FRAC_PI_2)),
            f64::EPSILON
        ));
        assert_eq!(2.0, gc2e_distance(Radians(core::f64::consts::PI)));
    }

    #[test]
    fn test_distance_and_azimuth_functions() {
        let angle_0 = Angle::default();
        let angle_45 = Angle::from_y_x(1.0, 1.0);

        let gc_distance = calculate_gc_distance(angle_0, angle_45, -angle_45);
        let haversine_distance =
            calculate_haversine_distance(angle_0, angle_45, -angle_45, -angle_45);
        assert!(is_within_tolerance(
            gc_distance.0,
            haversine_distance.0,
            f64::EPSILON
        ));

        let gc_azimuth = calculate_gc_azimuth(angle_0, angle_45, -angle_45);
        assert!(is_within_tolerance(
            -35.264389682754654,
            Degrees::from(gc_azimuth).0,
            32.0 * f64::EPSILON
        ));

        let gc_distance = calculate_gc_distance(angle_45, angle_45, angle_0);
        assert_eq!(0.0, gc_distance.0);
        let haversine_distance = calculate_haversine_distance(angle_45, angle_45, angle_0, angle_0);
        assert_eq!(0.0, haversine_distance.0);

        let gc_azimuth = calculate_gc_azimuth(angle_45, angle_45, angle_0);
        assert_eq!(0.0, Degrees::from(gc_azimuth).0);
    }

    #[test]
    fn test_north_and_south_pole_azimuths() {
        let angle_90 = Angle::new(trig::UnitNegRange(1.0), trig::UnitNegRange(0.0));
        let angle_m90 = -angle_90;

        let angle_0 = Angle::default();
        let angle_45 = Angle::from_y_x(1.0, 1.0);
        let angle_180 = Angle::new(trig::UnitNegRange(0.0), trig::UnitNegRange(-1.0));

        // From South Pole
        let gc_azimuth = calculate_gc_azimuth(angle_m90, angle_45, angle_0);
        assert_eq!(angle_0, gc_azimuth);

        // From North Pole
        let gc_azimuth = calculate_gc_azimuth(angle_90, angle_45, angle_0);
        assert_eq!(angle_180, gc_azimuth);
    }

    #[test]
    fn test_calculate_other_azimuth() {
        let angle_90 = Angle::from(Degrees(90.0));
        let angle_50 = Angle::from(Degrees(50.0));
        let angle_45 = Angle::from(Degrees(45.0));
        let angle_20 = Angle::from(Degrees(20.0));

        let result: Angle = calculate_other_azimuth(angle_20, angle_50, angle_20);
        assert!(is_within_tolerance(
            30.0,
            Degrees::from(result).0,
            16.0 * f64::EPSILON
        ));

        let result: Angle = calculate_other_azimuth(angle_50, angle_20, angle_20);
        assert!(is_within_tolerance(
            13.530064432438888,
            Degrees::from(result).0,
            16.0 * f64::EPSILON
        ));

        let result: Angle = calculate_other_azimuth(-angle_50, angle_50, angle_20);
        assert_eq!(20.0, Degrees::from(result).0);

        let result: Angle = calculate_other_azimuth(angle_45, angle_45, angle_90);
        assert_eq!(90.0, Degrees::from(result).0);
    }

    #[test]
    fn test_calculate_lat_delta_lon() {}
}
