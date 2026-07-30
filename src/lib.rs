// Copyright (c) 2024-2026 Via Technology Ltd.

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

//! [![crates.io](https://img.shields.io/crates/v/icao-isa.svg)](https://crates.io/crates/icao-isa)
//! [![docs.io](https://docs.rs/icao-isa/badge.svg)](https://docs.rs/icao-isa/)
//! [![License](https://img.shields.io/badge/License-MIT-blue)](https://opensource.org/license/mit/)
//! [![Rust](https://github.com/kenba/icao-isa-rs/actions/workflows/rust.yml/badge.svg)](https://github.com/kenba/icao-isa-rs/actions)
//! [![codecov](https://codecov.io/gh/kenba/icao-isa-rs/graph/badge.svg?token=6DTOY9Y4BT)](https://codecov.io/gh/kenba/icao-isa-rs)
//!
//! An implementation of the [International Civil Aviation Organization](https://icao.int/) (ICAO)
//! [International Standard Atmosphere](https://en.wikipedia.org/wiki/International_Standard_Atmosphere)
//! (ISA), see [ICAO Doc 7488/3](https://standart.aero/en/icao/book/doc-7488-manual-of-the-icao-standard-atmosphere-extended-to-80-kilometres-262-500-feet-en-cons).
//!
//! The library also includes functions for calculating:
//!
//! - true airspeed ([TAS](https://en.wikipedia.org/wiki/True_airspeed)) from calibrated airspeed ([CAS](https://en.wikipedia.org/wiki/Calibrated_airspeed)), pressure and temperature;
//! - CAS from TAS, pressure and temperature;
//! - TAS from [Mach number](https://en.wikipedia.org/wiki/Mach_number) and temperature;
//! - and the crossover altitude between CAS / MACH flight regimes.
//!
//! The equations for the functions above are from
//! [BADA User Manual revision 3-12](https://www.scribd.com/document/289480324/1-User-Manual-Bada-3-12).
//!
//! The library is declared [no_std](https://docs.rust-embedded.org/book/intro/no-std.html)
//! so it can be used in embedded applications.

#![cfg_attr(not(test), no_std)]
#![allow(clippy::suboptimal_flops)]

pub mod constants;

use icao_units::si::{
    Kelvin, KilogramsPerCubicMetre, Metres, MetresPerSecond, MetresPerSecondSquared, Pascals,
};
use num_traits::Float;

/// The pressure at `ISA_TROPOPAUSE_ALTITUDE` in `Pascals`.
/// See BADA Equation Eq 3.1-19
const ISA_TROPOPAUSE_PRESSURE: f64 = 22_632.040_095_007_81;

#[allow(clippy::missing_panics_doc)]
#[must_use]
pub fn g<T: Float>() -> MetresPerSecondSquared<T> {
    let value = T::from(constants::G).expect("Could not convert constant to Float");
    MetresPerSecondSquared(value)
}

#[allow(clippy::missing_panics_doc)]
#[must_use]
pub fn k<T: Float>() -> T {
    T::from(constants::K).expect("Could not convert constant to Float")
}

#[allow(clippy::missing_panics_doc)]
#[must_use]
pub fn r<T: Float>() -> T {
    T::from(constants::R).expect("Could not convert constant to Float")
}

#[allow(clippy::missing_panics_doc)]
#[must_use]
pub fn isa_sea_level_temperature<T: Float>() -> Kelvin<T> {
    let value =
        T::from(constants::ISA_SEA_LEVEL_TEMPERATURE).expect("Could not convert constant to Float");
    Kelvin(value)
}

#[allow(clippy::missing_panics_doc)]
#[must_use]
pub fn isa_sea_level_pressure<T: Float>() -> Pascals<T> {
    let value =
        T::from(constants::ISA_SEA_LEVEL_PRESSURE).expect("Could not convert constant to Float");
    Pascals(value)
}

#[allow(clippy::missing_panics_doc)]
#[must_use]
pub fn isa_sea_level_speed_of_sound<T: Float>() -> MetresPerSecond<T> {
    let value = T::from(constants::ISA_SEA_LEVEL_SPEED_OF_SOUND)
        .expect("Could not convert constant to Float");
    MetresPerSecond(value)
}

#[allow(clippy::missing_panics_doc)]
#[must_use]
pub fn isa_temperature_gradient<T: Float>() -> T {
    T::from(constants::ISA_TEMPERATURE_GRADIENT).expect("Could not convert constant to Float")
}

#[allow(clippy::missing_panics_doc)]
#[must_use]
pub fn isa_tropopause_temperature<T: Float>() -> Kelvin<T> {
    let value = T::from(constants::ISA_TROPOPAUSE_TEMPERATURE)
        .expect("Could not convert constant to Float");
    Kelvin(value)
}

#[allow(clippy::missing_panics_doc)]
#[must_use]
pub fn isa_tropopause_altitude<T: Float>() -> Metres<T> {
    let value =
        T::from(constants::ISA_TROPOPAUSE_ALTITUDE).expect("Could not convert constant to Float");
    Metres(value)
}

#[allow(clippy::missing_panics_doc)]
#[must_use]
pub fn isa_tropopause_pressure<T: Float>() -> Pascals<T> {
    let value = T::from(ISA_TROPOPAUSE_PRESSURE).expect("Could not convert constant to Float");
    Pascals(value)
}

/// The coefficient used in CAS / TAS conversions.
/// See BADA Equation 3.2-14
#[must_use]
fn u<T: Float>() -> T {
    const U: f64 = (constants::K - 1.0) / constants::K;
    T::from(U).expect("Could not convert constant to Float")
}

/// Another coefficient used in pressure conversions.
/// See BADA Equation 3.2-14
#[must_use]
fn inv_u<T: Float>() -> T {
    const INV_U: f64 = constants::K / (constants::K - 1.0);
    T::from(INV_U).expect("Could not convert constant to Float")
}

#[must_use]
fn k_minus_1_over_2<T: Float>() -> T {
    const K_MINUS_1_OVER_2: f64 = (constants::K - 1.0) / 2.0;
    T::from(K_MINUS_1_OVER_2).expect("Could not convert constant to Float")
}

#[must_use]
fn pressure_power<T: Float>() -> T {
    const PRESSURE_POWER: f64 =
        -constants::G / (constants::ISA_TEMPERATURE_GRADIENT * constants::R);
    T::from(PRESSURE_POWER).expect("Could not convert constant to Float")
}

#[must_use]
fn temperature_power<T: Float>() -> T {
    const TEMPERATURE_POWER: f64 =
        (constants::ISA_TEMPERATURE_GRADIENT * constants::R) / -constants::G;
    T::from(TEMPERATURE_POWER).expect("Could not convert constant to Float")
}

/// The factor used in calculating the density and pressure above Tropopause.
/// See BADA Equation 3.2-16
#[must_use]
fn tropopause_pressure_factor<T: Float>() -> T {
    const TROPOPAUSE_PRESSURE_FACTOR: f64 =
        -constants::G / (constants::R * constants::ISA_TROPOPAUSE_TEMPERATURE);
    T::from(TROPOPAUSE_PRESSURE_FACTOR).expect("Could not convert constant to Float")
}

/// Calculate the ISA pressure below the tropopause for the given altitude.
///
/// See BADA Rev 3.12, Eq 3.1-18
/// * `altitude` the altitude in Metres: max tropopause - `11_000` metres.
///
/// returns the pressure in pascals.
#[must_use]
fn calculate_troposphere_pressure<T: Float>(altitude: Metres<T>) -> Pascals<T> {
    const TEMPERATURE_FACTOR: f64 =
        constants::ISA_TEMPERATURE_GRADIENT / constants::ISA_SEA_LEVEL_TEMPERATURE;
    let temperature_factor =
        T::from(TEMPERATURE_FACTOR).expect("Could not convert constant to Float");

    let value = isa_sea_level_pressure::<T>().0
        * (T::one() + altitude.0 * temperature_factor).powf(pressure_power());
    Pascals(value)
}

/// Calculate the ISA pressure in the tropopause for the given altitude.
///
/// See BADA Rev 3.12, Eq 3.1-20
/// * `altitude` the altitude in Metres: min tropopause - `11_000` metres.
///
/// returns the pressure in Pascals.
#[must_use]
fn calculate_tropopause_pressure<T: Float>(altitude: Metres<T>) -> Pascals<T> {
    let value = isa_tropopause_pressure::<T>().0
        * (tropopause_pressure_factor::<T>() * (altitude.0 - isa_tropopause_altitude().0)).exp();
    Pascals(value)
}

/// Calculate the ISA pressure corresponding to the given altitude.
///
/// Note: ISA pressure does **NOT** vary with temperature.
/// See BADA Rev 3.12, Eq 3.1-18 & Eq 3.1-20
/// * `altitude` the pressure altitude in metres.
///
/// returns the pressure in Pascals.
#[must_use]
pub fn calculate_isa_pressure<T: Float>(altitude: Metres<T>) -> Pascals<T> {
    if altitude < isa_tropopause_altitude() {
        calculate_troposphere_pressure(altitude)
    } else {
        calculate_tropopause_pressure(altitude)
    }
}

/// Calculate the altitude corresponding to the given pressure below the tropopause.
///
/// See BADA Rev 3.12, Eq 3.1-18
/// * `pressure` the pressure in Pascals.
///
/// returns the altitude in metres.
#[must_use]
fn calculate_troposphere_altitude<T: Float>(pressure: Pascals<T>) -> Metres<T> {
    let pressure_ratio = pressure.0 / isa_sea_level_pressure().0;
    let altitude_ratio = pressure_ratio.powf(temperature_power()) - T::one();
    let value = altitude_ratio * isa_sea_level_temperature().0 / isa_temperature_gradient();
    Metres(value)
}

/// Calculate the altitude corresponding to the given pressure in the tropopause.
///
/// See BADA Rev 3.12, Eq 3.1-20
/// * `pressure` the pressure in Pascals.
///
/// returns the altitude in metres.
#[must_use]
fn calculate_tropopause_altitude<T: Float>(pressure: Pascals<T>) -> Metres<T> {
    let altitude_delta =
        (pressure.0 / isa_tropopause_pressure().0).ln() / tropopause_pressure_factor();
    let value = altitude_delta + isa_tropopause_altitude().0;
    Metres(value)
}

/// Calculate the ISA altitude corresponding to the given pressure.
///
/// See BADA Rev 3.12, Eq 3.1-18 & Eq 3.1-20
/// * `altitude` the pressure altitude in metres.
///
/// returns the pressure in Pascals.
#[must_use]
pub fn calculate_isa_altitude<T: Float>(pressure: Pascals<T>) -> Metres<T> {
    if pressure > isa_tropopause_pressure() {
        calculate_troposphere_altitude(pressure)
    } else {
        calculate_tropopause_altitude(pressure)
    }
}

/// Calculate the ISA temperature corresponding to the given altitude and
/// difference in Sea level temperature.
///
/// See ICAO Doc 7488/3, Eq (11)
/// * `altitude` the altitude in Metres.
/// * `delta_temperature` the difference from ISA temperature at Sea level.
///
/// returns the temperature in Kelvin.
#[must_use]
pub fn calculate_isa_temperature<T: Float>(
    altitude: Metres<T>,
    delta_temperature: Kelvin<T>,
) -> Kelvin<T> {
    let temperature = Kelvin(
        delta_temperature.0
            + isa_sea_level_temperature().0
            + altitude.0 * isa_temperature_gradient(),
    );

    if temperature > isa_tropopause_temperature() {
        temperature
    } else {
        isa_tropopause_temperature()
    }
}

/// Calculate the air density given the air temperature and pressure.\
/// Uses the Ideal Gas Equation (Boyles law)
///
/// See See ICAO Doc 7488/3, Eq (3).
/// * `pressure` the pressure in Pascals.
/// * `temperature` the temperature in Kelvin.
///
/// returns the density in Kg per cubic metre.
#[must_use]
pub fn calculate_density<T: Float>(
    pressure: Pascals<T>,
    temperature: Kelvin<T>,
) -> KilogramsPerCubicMetre<T> {
    KilogramsPerCubicMetre(pressure.0 / (temperature.0 * r()))
}

/// Calculate the True Air Speed (TAS) from the Calibrated Air Speed (CAS)
/// at the given pressure and temperature.
///
/// See BADA Rev 3.12, Eq 3.1.23
/// * `cas` the Calibrated Air Speed in metres per second.
/// * `pressure` the pressure in Pascals.
/// * `temperature` the temperature in Kelvin.
///
/// returns the True Air Speed in metres per second.
#[must_use]
pub fn calculate_true_air_speed<T: Float>(
    cas: MetresPerSecond<T>,
    pressure: Pascals<T>,
    temperature: Kelvin<T>,
) -> MetresPerSecond<T> {
    let two_r = r::<T>() + r::<T>();
    let outer_factor = two_r / u();
    let inner_factor = T::one() / (outer_factor * isa_sea_level_temperature().0);

    let cas_factor = (T::one() + inner_factor * cas.0 * cas.0).powf(inv_u()) - T::one();
    let cas_pressure_factor =
        (T::one() + cas_factor * isa_sea_level_pressure().0 / pressure.0).powf(u()) - T::one();
    let value = (outer_factor * temperature.0 * cas_pressure_factor).sqrt();
    MetresPerSecond(value)
}

/// Calculate the Calibrated Air Speed (CAS) from the True Air Speed (TAS)
/// at the given pressure and temperature.
///
/// See BADA Rev 3.12, Eq 3.1.24
/// * `tas` the True Air Speed in metres per second.
/// * `pressure` the pressure in Pascals.
/// * `temperature` the temperature in Kelvin.
///
/// * returns the Calibrated Air Speed in metres per second.
#[must_use]
pub fn calculate_calibrated_air_speed<T: Float>(
    tas: MetresPerSecond<T>,
    pressure: Pascals<T>,
    temperature: Kelvin<T>,
) -> MetresPerSecond<T> {
    let two_r = r::<T>() + r::<T>();
    let inner_factor = u::<T>() / two_r;
    let outer_factor = isa_sea_level_temperature::<T>().0 / inner_factor;

    let tas_factor =
        (T::one() + inner_factor * tas.0 * tas.0 / temperature.0).powf(inv_u()) - T::one();
    let tas_pressure_factor =
        (T::one() + pressure.0 * tas_factor / isa_sea_level_pressure().0).powf(u()) - T::one();
    let value = (outer_factor * tas_pressure_factor).sqrt();
    MetresPerSecond(value)
}

/// Calculate the speed of sound for the given temperature.
///
/// See ICAO Doc 7488/3, Eq 21
/// * `temperature` the temperature in Kelvin.
///
/// returns the speed of sound in metres per second.
#[allow(clippy::missing_panics_doc)]
#[must_use]
pub fn speed_of_sound<T: Float>(temperature: Kelvin<T>) -> MetresPerSecond<T> {
    const KR: f64 = constants::K * constants::R;
    let kr = T::from(KR).expect("Could not convert constant to Float");

    MetresPerSecond((temperature.0 * kr).sqrt())
}

/// Calculate the True Air Speed (TAS) from the Mach number at the given temperature.
///
/// See BADA Rev 3.12, Eq 3.1.22
/// * `mach` the Mach number.
/// * `temperature` the temperature in Kelvin.
///
/// returns the True Air Speed in metres per second.
#[must_use]
pub fn mach_true_air_speed<T: Float>(mach: T, temperature: Kelvin<T>) -> MetresPerSecond<T> {
    MetresPerSecond(mach * speed_of_sound(temperature).0)
}

/// This function calculates the crossover pressure ratio between the
/// Calibrated Air Speed (CAS) and Mach number.
///
/// See BADA Rev 3.12, Eq 3.1-29
/// * `cas` the Calibrated Air Speed in metres per second.
/// * `mach` the Mach number.
///
/// returns the crossover pressure ratio.
#[must_use]
fn calculate_crossover_pressure_ratio<T: Float>(cas: MetresPerSecond<T>, mach: T) -> T {
    let cas_mach = cas.0 / isa_sea_level_speed_of_sound().0;
    let numerator =
        (T::one() + k_minus_1_over_2::<T>() * cas_mach * cas_mach).powf(inv_u()) - T::one();
    let denominator = (T::one() + k_minus_1_over_2::<T>() * mach * mach).powf(inv_u()) - T::one();

    numerator / denominator
}

/// Calculate the crossover altitude at which the True Air Speeds (TAS)
/// corresponding to the given Calibrated Air Speed (CAS) and Mach number are
/// the same.
///
/// See BADA Rev 3.12, Eq 3.1-27
/// * `cas` the Calibrated Air Speed in metres per second.
/// * `mach` the Mach number.
///
/// returns the altitude in metres.
#[must_use]
pub fn calculate_crossover_altitude<T: Float>(cas: MetresPerSecond<T>, mach: T) -> Metres<T> {
    let isa_temperature_ratio =
        isa_sea_level_temperature::<T>().0 / -isa_temperature_gradient::<T>();

    let temperature_ratio = calculate_crossover_pressure_ratio(cas, mach).powf(temperature_power());
    let value = isa_temperature_ratio * (T::one() - temperature_ratio);
    Metres(value)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_constant_functions() {
        assert_eq!(MetresPerSecondSquared(constants::G), g());
        assert_eq!(constants::K, k());
        assert_eq!(constants::R, r());
    }

    #[test]
    fn test_calculate_isa_pressure() {
        // calculate_troposphere_pressure
        assert_eq!(
            Pascals(constants::ISA_SEA_LEVEL_PRESSURE),
            calculate_isa_pressure(Metres(0.0))
        );
        assert!((89874.563 - calculate_isa_pressure(Metres(1000.0)).0).abs() < 0.001);
        assert!((79495.201 - calculate_isa_pressure(Metres(2000.0)).0).abs() < 0.001);
        assert!((22635.609 - calculate_isa_pressure(Metres(10999.0)).0).abs() < 0.001);

        // calculate_tropopause_pressure
        assert_eq!(
            22632.04009500781,
            calculate_isa_pressure(Metres(constants::ISA_TROPOPAUSE_ALTITUDE)).0
        );
        assert!((19330.383 - calculate_isa_pressure(Metres(12000.0)).0).abs() < 0.001);
    }

    #[test]
    fn test_calculate_isa_altitude() {
        // calculate_troposphere_altitude
        assert_eq!(
            0.0,
            calculate_isa_altitude(Pascals(constants::ISA_SEA_LEVEL_PRESSURE)).0
        );
        assert!((1000.0 - calculate_isa_altitude(Pascals(89874.563)).0).abs() < 0.001);
        assert!((2000.0 - calculate_isa_altitude(Pascals(79495.201)).0).abs() < 0.001);
        assert!((10999.0 - calculate_isa_altitude(Pascals(22635.609)).0).abs() < 0.001);

        // calculate_tropopause_altitude
        assert_eq!(
            constants::ISA_TROPOPAUSE_ALTITUDE,
            calculate_isa_altitude(Pascals(ISA_TROPOPAUSE_PRESSURE)).0
        );
        assert!((12000.0 - calculate_isa_altitude(Pascals(19330.383)).0).abs() < 0.001);
    }

    #[test]
    fn test_calculate_isa_temperature() {
        assert_eq!(
            constants::ISA_SEA_LEVEL_TEMPERATURE - 3.25,
            calculate_isa_temperature(Metres(500.0), Kelvin(0.0)).0
        );
        assert_eq!(
            constants::ISA_SEA_LEVEL_TEMPERATURE - 13.0,
            calculate_isa_temperature(Metres(2000.0), Kelvin(0.0)).0
        );
        assert_eq!(
            constants::ISA_TROPOPAUSE_TEMPERATURE,
            calculate_isa_temperature(Metres(constants::ISA_TROPOPAUSE_ALTITUDE), Kelvin(0.0)).0
        );
        assert!(
            (constants::ISA_TROPOPAUSE_TEMPERATURE + 10.
                - calculate_isa_temperature(
                    Metres(constants::ISA_TROPOPAUSE_ALTITUDE),
                    Kelvin(10.0)
                )
                .0)
                .abs()
                < 1.0e-9
        );
        assert_eq!(
            constants::ISA_TROPOPAUSE_TEMPERATURE,
            calculate_isa_temperature(Metres(12000.0), Kelvin(-10.0)).0
        );
    }

    #[test]
    fn test_calculate_density() {
        assert!(
            (constants::ISA_SEA_LEVEL_DENSITY
                - calculate_density(
                    Pascals(constants::ISA_SEA_LEVEL_PRESSURE),
                    Kelvin(constants::ISA_SEA_LEVEL_TEMPERATURE)
                )
                .0)
                .abs()
                < 2.0e-8
        );
        assert!(
            (0.3639176
                - calculate_density(
                    Pascals(ISA_TROPOPAUSE_PRESSURE),
                    Kelvin(constants::ISA_TROPOPAUSE_TEMPERATURE)
                )
                .0)
                .abs()
                < 1.0e-6
        );
    }

    #[test]
    fn test_calculate_true_air_speed() {
        assert!(
            (150.0
                - calculate_true_air_speed(
                    MetresPerSecond(150.0),
                    Pascals(constants::ISA_SEA_LEVEL_PRESSURE),
                    Kelvin(constants::ISA_SEA_LEVEL_TEMPERATURE)
                )
                .0)
                .abs()
                < 1.0e-9
        );
        assert!(
            (164.458
                - calculate_true_air_speed(
                    MetresPerSecond(150.0),
                    Pascals(79495.201),
                    Kelvin(constants::ISA_SEA_LEVEL_TEMPERATURE - 13.0)
                )
                .0)
                .abs()
                < 0.001
        );
    }

    #[test]
    fn test_calculate_calibrated_air_speed() {
        assert!(
            (150.0
                - calculate_calibrated_air_speed(
                    MetresPerSecond(150.0),
                    Pascals(constants::ISA_SEA_LEVEL_PRESSURE),
                    Kelvin(constants::ISA_SEA_LEVEL_TEMPERATURE)
                )
                .0)
                .abs()
                < 1.0e-9
        );
        assert!(
            (150.0
                - calculate_calibrated_air_speed(
                    MetresPerSecond(164.458),
                    Pascals(79495.201),
                    Kelvin(constants::ISA_SEA_LEVEL_TEMPERATURE - 13.0)
                )
                .0)
                .abs()
                < 0.001
        );
    }

    #[test]
    fn test_speed_of_sound() {
        assert_eq!(0.0, speed_of_sound(Kelvin(0.0)).0);
        assert!(
            (constants::ISA_SEA_LEVEL_SPEED_OF_SOUND
                - speed_of_sound(Kelvin(constants::ISA_SEA_LEVEL_TEMPERATURE)).0)
                .abs()
                < 0.001
        );
        assert!(
            (295.070 - speed_of_sound(Kelvin(constants::ISA_TROPOPAUSE_TEMPERATURE)).0).abs()
                < 0.001
        );
    }

    #[test]
    fn test_mach_true_air_speed() {
        assert!(
            (0.8 * constants::ISA_SEA_LEVEL_SPEED_OF_SOUND
                - mach_true_air_speed(0.8, Kelvin(constants::ISA_SEA_LEVEL_TEMPERATURE)).0)
                .abs()
                < 0.001
        );
        assert!(
            (250.809 - mach_true_air_speed(0.85, Kelvin(constants::ISA_TROPOPAUSE_TEMPERATURE)).0)
                .abs()
                < 0.001
        );
    }

    #[test]
    fn test_calculate_crossover_altitude() {
        let cas = MetresPerSecond(155.0);
        let crossover_altitude = calculate_crossover_altitude(cas, 0.79);
        assert!((9070.814 - crossover_altitude.0).abs() < 0.001);

        // The TAS should be the same from both CAS and MACH at the crossover_altitude
        let _pressure = calculate_isa_pressure(crossover_altitude);
        let _temperature = calculate_isa_temperature(crossover_altitude, Kelvin(0.0));
        // let tas_from_cas = calculate_true_air_speed(cas, pressure, temperature);
        // let tas_from_mach = mach_true_air_speed(0.79, temperature);
        // assert!((tas_from_cas.0 - tas_from_mach.0).abs() < 0.001);
    }
}
