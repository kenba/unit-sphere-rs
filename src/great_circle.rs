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

//! The `great_circle` module contains functions for calculating the course
//! and distance between points along great circles on a unit sphere.

use angle_sc::{Angle, Radians, trig};
use num_traits::{Float, float::FloatConst};

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
pub fn calculate_haversine_distance<T: Float + FloatConst>(
    a_lat: Angle<T>,
    b_lat: Angle<T>,
    delta_long: Angle<T>,
    delta_lat: Angle<T>,
) -> Radians<T> {
    let min_value = T::epsilon() + T::epsilon();
    let haversine_lat = trig::sq_sine_half(delta_lat.cos());
    let haversine_lon = trig::sq_sine_half(delta_long.cos());

    let a =
        (haversine_lat + a_lat.cos().0 * b_lat.cos().0 * haversine_lon).clamp(T::zero(), T::one());

    if a < min_value {
        Radians::<T>(T::zero())
    } else {
        let half_value = a.sqrt().asin();
        Radians::<T>(half_value + half_value)
    }
}

/// Convert a Euclidean distance to a Great Circle distance (in `Radians`).
/// e should satisfy: 0 <= e <= 2, if not it is clamped into range.
#[must_use]
pub fn e2gc_distance<T: Float + FloatConst>(e: T) -> Radians<T> {
    let two = T::one() + T::one();
    let half = T::one() / two;
    Radians::<T>(two * (half * e.clamp(T::zero(), two)).asin())
}

/// Convert a Great Circle distance (in radians) to a Euclidean distance.
#[must_use]
pub fn gc2e_distance<T: Float + FloatConst>(gc: Radians<T>) -> T {
    let two = T::one() + T::one();
    let half = T::one() / two;
    two * (half * gc.0).sin()
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
pub fn sq_euclidean_distance<T: Float>(
    a_lat: Angle<T>,
    b_lat: Angle<T>,
    delta_long: Angle<T>,
) -> T {
    let two = T::one() + T::one();
    let four = two + two;

    let delta_x = b_lat.cos().0 * delta_long.cos().0 - a_lat.cos().0;
    let delta_y = b_lat.cos().0 * delta_long.sin().0;
    let delta_z = b_lat.sin().0 - a_lat.sin().0;

    let result = delta_x * delta_x + delta_y * delta_y + delta_z * delta_z;
    result.clamp(T::zero(), four)
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
pub fn calculate_gc_distance<T: Float + FloatConst>(
    a_lat: Angle<T>,
    b_lat: Angle<T>,
    delta_long: Angle<T>,
) -> Radians<T> {
    e2gc_distance(sq_euclidean_distance(a_lat, b_lat, delta_long).sqrt())
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
pub fn calculate_gc_azimuth<T: Float>(
    a_lat: Angle<T>,
    b_lat: Angle<T>,
    delta_long: Angle<T>,
) -> Angle<T> {
    let min_value = T::epsilon() + T::epsilon();

    // if start point is North or South pole
    if a_lat.cos().0 < min_value {
        if a_lat.sin().0.is_sign_negative() {
            // South pole, azimuth is zero
            Angle::default()
        } else {
            // North pole, azimuth is 180 degrees
            Angle::new(trig::UnitNegRange(T::zero()), trig::UnitNegRange(-T::one()))
        }
    } else {
        let sin_azimuth = b_lat.cos().0 * delta_long.sin().0;
        let temp = (a_lat.sin().0 * b_lat.cos().0 * delta_long.sin().0 * delta_long.sin().0)
            / (T::one() + delta_long.cos().abs().0);
        let cos_azimuth = if delta_long.cos().0 < T::zero() {
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
#[must_use]
pub fn calculate_sigma<T: Float>(
    a_lat: Angle<T>,
    b_lat: Angle<T>,
    delta_long: Angle<T>,
) -> Angle<T> {
    let a = b_lat.cos().0 * delta_long.sin().0;
    let b = a_lat.cos().0 * b_lat.sin().0 - a_lat.sin().0 * b_lat.cos().0 * delta_long.cos().0;
    let sin_sigma = trig::UnitNegRange::<T>(a.hypot(b));
    let cos_sigma = trig::UnitNegRange::<T>(
        a_lat.sin().0 * b_lat.sin().0 + a_lat.cos().0 * b_lat.cos().0 * delta_long.cos().0,
    );
    Angle::new(sin_sigma, cos_sigma)
}

/// Calculate the latitude at great circle distance, sigma.
///
/// From Eurocae ED-323 Appendix E.
/// * `a_lat` - the latitude of the start point.
/// * `azi` - the azimuth at `a_lat`.
/// * `sigma` - the distance on the auxiliary sphere as an Angle.
///
/// return the latitude at sigma from `a_lat`.
#[must_use]
pub fn calculate_latitude<T: Float>(lat: Angle<T>, azi: Angle<T>, sigma: Angle<T>) -> Angle<T> {
    let sin_lat = trig::UnitNegRange::<T>(
        lat.sin().0 * sigma.cos().0 + lat.cos().0 * sigma.sin().0 * azi.cos().0,
    );
    let cos_lat = trig::UnitNegRange::<T>(
        (lat.cos().0 * azi.sin().0)
            .hypot(lat.sin().0 * sigma.sin().0 - lat.cos().0 * sigma.cos().0 * azi.cos().0),
    );

    Angle::new(sin_lat, cos_lat)
}

/// Calculate the longitude difference at great circle distance, sigma.
/// * `a_lat` - the latitude of the start point.
/// * `azi` - the azimuth at `a_lat`.
/// * `sigma` - the distance on the auxiliary sphere as an Angle.
///
/// return the longitude difference from `a_lat`.
#[must_use]
pub fn calculate_delta_longitude<T: Float>(
    lat: Angle<T>,
    azi: Angle<T>,
    sigma: Angle<T>,
) -> Angle<T> {
    let sin_lon = sigma.sin().0 * azi.sin().0;
    let cos_lon = lat.cos().0 * sigma.cos().0 - lat.sin().0 * sigma.sin().0 * azi.cos().0;

    Angle::from_y_x(sin_lon, cos_lon)
}

#[cfg(test)]
mod tests {
    use super::*;
    use angle_sc::{Degrees, is_within_tolerance};

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
        let sigma = calculate_sigma(angle_0, angle_45, -angle_45);
        assert!(is_within_tolerance(
            gc_distance.0,
            Radians::from(sigma).0,
            f64::EPSILON
        ));

        let gc_azimuth = calculate_gc_azimuth(angle_0, angle_45, -angle_45);
        assert!(is_within_tolerance(
            -35.264389682754654,
            Degrees::from(gc_azimuth).0,
            32.0 * f64::EPSILON
        ));

        let latitude = calculate_latitude(angle_0, gc_azimuth, sigma);
        assert_eq!(angle_45, latitude);

        let delta_lon = calculate_delta_longitude(angle_0, gc_azimuth, sigma);
        assert!(is_within_tolerance(
            -45.0,
            Degrees::from(delta_lon).0,
            f64::EPSILON
        ));

        let gc_distance = calculate_gc_distance(angle_45, angle_45, angle_0);
        assert_eq!(0.0, gc_distance.0);
        let haversine_distance = calculate_haversine_distance(angle_45, angle_45, angle_0, angle_0);
        assert_eq!(0.0, haversine_distance.0);
        let sigma = calculate_sigma(angle_45, angle_45, angle_0);
        assert_eq!(0.0, Radians::from(sigma).0);

        let gc_azimuth = calculate_gc_azimuth(angle_45, angle_45, angle_0);
        assert_eq!(0.0, Degrees::from(gc_azimuth).0);

        let latitude = calculate_latitude(angle_45, gc_azimuth, sigma);
        assert!(is_within_tolerance(
            45.0,
            Degrees::from(latitude).0,
            f64::EPSILON
        ));

        let delta_lon = calculate_delta_longitude(angle_45, gc_azimuth, sigma);
        assert_eq!(0.0, Degrees::from(delta_lon).0);
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
}
