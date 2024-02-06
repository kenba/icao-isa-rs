// Copyright (c) 2024 Via Technology Ltd.

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

use icao_units::si::{Kelvin, KilogramsPerCubicMetre, Metres, MetresPerSecond, Pascals};

/// The coefficient used in CAS / TAS conversions.
/// See BADA Equation 3.2-14
const U: f64 = (constants::K - 1.0) / constants::K;

/// Another coefficient used in pressure conversions.
/// See BADA Equation 3.2-14
const INV_U: f64 = 1.0 / U;

const K_MINUS_1_OVER_2: f64 = (constants::K - 1.0) / 2.0;

/// The Power factor used in calculating the pressure below the Tropopause.
/// See BADA Equation 3.1-18
const PRESSURE_POWER: f64 = -constants::g.0 / (constants::TEMPERATURE_GRADIENT * constants::R);

/// The Power factor used in calculating the altitude below the Tropopause.
/// See BADA Equation 3.1-8
const TEMPERATURE_POWER: f64 = 1.0 / PRESSURE_POWER;

/// The pressure at `TROPOPAUSE_ALTITUDE`.
/// See BADA Equation Eq 3.1-19
const TROPOPAUSE_PRESSURE: Pascals = Pascals(22_632.040_095_007_81);

/// The factor used in calculating the density and pressure above Tropopause.
/// See BADA Equation 3.2-16
const TROPOPAUSE_PRESSURE_FACTOR: f64 =
    -constants::g.0 / (constants::R * constants::TROPOPAUSE_TEMPERATURE.0);

/// Calculate the ISA pressure below the tropopause for the given altitude.  
/// See BADA Rev 3.12, Eq 3.1-18
/// * `altitude` the altitude in Metres: max tropopause - `11_000` metres.
///
/// returns the pressure in pascals.
#[must_use]
fn calculate_troposphere_pressure(altitude: Metres) -> Pascals {
    Pascals(
        constants::SEA_LEVEL_PRESSURE.0
            * libm::pow(
                1.0 + altitude.0 * constants::TEMPERATURE_GRADIENT
                    / constants::SEA_LEVEL_TEMPERATURE.0,
                PRESSURE_POWER,
            ),
    )
}

/// Calculate the ISA pressure in the tropopause for the given altitude.  
/// See BADA Rev 3.12, Eq 3.1-20
/// * `altitude` the altitude in Metres: min tropopause - `11_000` metres.
///
/// returns the pressure in Pascals.
#[must_use]
fn calculate_tropopause_pressure(altitude: Metres) -> Pascals {
    Pascals(
        TROPOPAUSE_PRESSURE.0
            * libm::exp(
                TROPOPAUSE_PRESSURE_FACTOR * (altitude.0 - constants::TROPOPAUSE_ALTITUDE.0),
            ),
    )
}

/// Calculate the ISA pressure corresponding to the given altitude.  
/// Note: ISA pressure does **NOT** vary with temperature.  
/// See BADA Rev 3.12, Eq 3.1-18 & Eq 3.1-20
/// * `altitude` the pressure altitude in metres.
///
/// returns the pressure in Pascals.
#[must_use]
pub fn calculate_isa_pressure(altitude: Metres) -> Pascals {
    if altitude < constants::TROPOPAUSE_ALTITUDE {
        calculate_troposphere_pressure(altitude)
    } else {
        calculate_tropopause_pressure(altitude)
    }
}

/// Calculate the altitude corresponding to the given pressure below the tropopause.  
/// See BADA Rev 3.12, Eq 3.1-18
/// * `pressure` the pressure in Pascals.
///
/// returns the altitude in metres.
#[must_use]
fn calculate_troposphere_altitude(pressure: Pascals) -> Metres {
    let pressure_ratio = pressure.0 / constants::SEA_LEVEL_PRESSURE.0;
    let altitude_ratio = libm::pow(pressure_ratio, TEMPERATURE_POWER) - 1.0;
    Metres(altitude_ratio * constants::SEA_LEVEL_TEMPERATURE.0 / constants::TEMPERATURE_GRADIENT)
}

/// Calculate the altitude corresponding to the given pressure in the tropopause.  
/// See BADA Rev 3.12, Eq 3.1-20
/// * `pressure` the pressure in Pascals.
///
/// returns the altitude in metres.
#[must_use]
fn calculate_tropopause_altitude(pressure: Pascals) -> Metres {
    let altitude_delta = libm::log(pressure.0 / TROPOPAUSE_PRESSURE.0) / TROPOPAUSE_PRESSURE_FACTOR;
    Metres(constants::TROPOPAUSE_ALTITUDE.0 + altitude_delta)
}

/// Calculate the ISA altitude corresponding to the given pressure.  
/// See BADA Rev 3.12, Eq 3.1-18 & Eq 3.1-20
/// * `altitude` the pressure altitude in metres.
///
/// returns the pressure in Pascals.
#[must_use]
pub fn calculate_isa_altitude(pressure: Pascals) -> Metres {
    if pressure > TROPOPAUSE_PRESSURE {
        calculate_troposphere_altitude(pressure)
    } else {
        calculate_tropopause_altitude(pressure)
    }
}

/// Calculate the ISA temperature corresponding to the given altitude and
/// difference in Sea level temperature.  
/// See ICAO Doc 7488/3, Eq (11)
/// * `altitude` the altitude in Metres.
/// * `delta_temperature` the difference from ISA temperature at Sea level.
///
/// returns the temperature in Kelvin.
#[must_use]
pub fn calculate_isa_temperature(altitude: Metres, delta_temperature: Kelvin) -> Kelvin {
    let temperature = Kelvin(
        constants::SEA_LEVEL_TEMPERATURE.0
            + delta_temperature.0
            + altitude.0 * constants::TEMPERATURE_GRADIENT,
    );

    if temperature > constants::TROPOPAUSE_TEMPERATURE {
        temperature
    } else {
        constants::TROPOPAUSE_TEMPERATURE
    }
}

/// Calculate the air density given the air temperature and pressure.  
/// Uses the Ideal Gas Equation (Boyles law)  
/// See See ICAO Doc 7488/3, Eq (3).
/// * `pressure` the pressure in Pascals.
/// * `temperature` the temperature in Kelvin.
///
/// returns the density in Kg per cubic metre.
#[must_use]
pub fn calculate_density(pressure: Pascals, temperature: Kelvin) -> KilogramsPerCubicMetre {
    KilogramsPerCubicMetre(pressure.0 / (temperature.0 * constants::R))
}

/// Calculate the True Air Speed (TAS) from the Calibrated Air Speed (CAS)
/// at the given pressure and temperature.  
/// See BADA Rev 3.12, Eq 3.1.23
/// * `cas` the Calibrated Air Speed in metres per second.
/// * `pressure` the pressure in Pascals.
/// * `temperature` the temperature in Kelvin.
///
/// returns the True Air Speed in metres per second.
#[must_use]
pub fn calculate_true_air_speed(
    cas: MetresPerSecond,
    pressure: Pascals,
    temperature: Kelvin,
) -> MetresPerSecond {
    const INNER_FACTOR: f64 = U / (2.0 * constants::R * constants::SEA_LEVEL_TEMPERATURE.0);
    const OUTER_FACTOR: f64 = 2.0 * constants::R / U;

    let cas_factor = libm::pow(1.0 + INNER_FACTOR * cas.0 * cas.0, INV_U) - 1.0;
    let cas_pressure_factor = libm::pow(
        1.0 + constants::SEA_LEVEL_PRESSURE.0 * cas_factor / pressure.0,
        U,
    ) - 1.0;
    MetresPerSecond(libm::sqrt(
        OUTER_FACTOR * temperature.0 * cas_pressure_factor,
    ))
}

/// Calculate the Calibrated Air Speed (CAS) from the True Air Speed (TAS)
/// at the given pressure and temperature.  
/// See BADA Rev 3.12, Eq 3.1.24
/// * `tas` the True Air Speed in metres per second.
/// * `pressure` the pressure in Pascals.
/// * `temperature` the temperature in Kelvin.
///
/// * returns the Calibrated Air Speed in metres per second.
#[must_use]
pub fn calculate_calibrated_air_speed(
    tas: MetresPerSecond,
    pressure: Pascals,
    temperature: Kelvin,
) -> MetresPerSecond {
    const INNER_FACTOR: f64 = U / (2.0 * constants::R);
    const OUTER_FACTOR: f64 = 2.0 * constants::R * constants::SEA_LEVEL_TEMPERATURE.0 / U;

    let tas_factor = libm::pow(1.0 + INNER_FACTOR * tas.0 * tas.0 / temperature.0, INV_U) - 1.0;
    let tas_pressure_factor = libm::pow(
        1.0 + pressure.0 * tas_factor / constants::SEA_LEVEL_PRESSURE.0,
        U,
    ) - 1.0;

    MetresPerSecond(libm::sqrt(OUTER_FACTOR * tas_pressure_factor))
}

/// Calculate the speed of sound for the given temperature.  
/// See ICAO Doc 7488/3, Eq 21
/// * `temperature` the temperature in Kelvin.
///
/// returns the speed of sound in metres per second.
#[must_use]
pub fn speed_of_sound(temperature: Kelvin) -> MetresPerSecond {
    MetresPerSecond(libm::sqrt(constants::K * constants::R * temperature.0))
}

/// Calculate the True Air Speed (TAS) from the Mach number at the given temperature.  
/// See BADA Rev 3.12, Eq 3.1.22
/// * `mach` the Mach number.
/// * `temperature` the temperature in Kelvin.
///
/// returns the True Air Speed in metres per second.
#[must_use]
pub fn mach_true_air_speed(mach: f64, temperature: Kelvin) -> MetresPerSecond {
    MetresPerSecond(mach * speed_of_sound(temperature).0)
}

/// This function calculates the crossover pressure ratio between the
/// Calibrated Air Speed (CAS) and Mach number.  
/// See BADA Rev 3.12, Eq 3.1-29
/// * `cas` the Calibrated Air Speed in metres per second.
/// * `mach` the Mach number.
///
/// returns the crossover pressure ratio.
#[must_use]
fn calculate_crossover_pressure_ratio(cas: MetresPerSecond, mach: f64) -> f64 {
    let cas_mach = cas.0 / constants::SEA_LEVEL_SPEED_OF_SOUND.0;

    let numerator = libm::pow(1.0 + K_MINUS_1_OVER_2 * cas_mach * cas_mach, INV_U) - 1.0;
    let denominator = libm::pow(1.0 + K_MINUS_1_OVER_2 * mach * mach, INV_U) - 1.0;

    numerator / denominator
}

/// Calculate the crossover altitude at which the True Air Speeds (TAS)
/// corresponding to the given Calibrated Air Speed (CAS) and Mach number are
/// the same.  
/// See BADA Rev 3.12, Eq 3.1-27
/// * `cas` the Calibrated Air Speed in metres per second.
/// * `mach` the Mach number.
///
/// returns the altitude in metres.
#[must_use]
pub fn calculate_crossover_altitude(cas: MetresPerSecond, mach: f64) -> Metres {
    let temperature_ratio = libm::pow(
        calculate_crossover_pressure_ratio(cas, mach),
        TEMPERATURE_POWER,
    );

    Metres(
        constants::SEA_LEVEL_TEMPERATURE.0 * (1.0 - temperature_ratio)
            / -constants::TEMPERATURE_GRADIENT,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_calculate_isa_pressure() {
        // calculate_troposphere_pressure
        assert_eq!(
            constants::SEA_LEVEL_PRESSURE,
            calculate_isa_pressure(Metres(0.0))
        );
        assert!(libm::fabs(89874.563 - calculate_isa_pressure(Metres(1000.0)).0) < 0.001);
        assert!(libm::fabs(79495.201 - calculate_isa_pressure(Metres(2000.0)).0) < 0.001);
        assert!(libm::fabs(22635.609 - calculate_isa_pressure(Metres(10999.0)).0) < 0.001);

        // calculate_tropopause_pressure
        assert_eq!(
            22632.04009500781,
            calculate_isa_pressure(constants::TROPOPAUSE_ALTITUDE).0
        );
        assert!(libm::fabs(19330.383 - calculate_isa_pressure(Metres(12000.0)).0) < 0.001);
    }

    #[test]
    fn test_calculate_isa_altitude() {
        // calculate_troposphere_altitude
        assert_eq!(0.0, calculate_isa_altitude(constants::SEA_LEVEL_PRESSURE).0);
        assert!(libm::fabs(1000.0 - calculate_isa_altitude(Pascals(89874.563)).0) < 0.001);
        assert!(libm::fabs(2000.0 - calculate_isa_altitude(Pascals(79495.201)).0) < 0.001);
        assert!(libm::fabs(10999.0 - calculate_isa_altitude(Pascals(22635.609)).0) < 0.001);

        // calculate_tropopause_altitude
        assert_eq!(
            constants::TROPOPAUSE_ALTITUDE.0,
            calculate_isa_altitude(TROPOPAUSE_PRESSURE).0
        );
        assert!(libm::fabs(12000.0 - calculate_isa_altitude(Pascals(19330.383)).0) < 0.001);
    }

    #[test]
    fn test_calculate_isa_temperature() {
        assert_eq!(
            constants::SEA_LEVEL_TEMPERATURE.0 - 3.25,
            calculate_isa_temperature(Metres(500.0), Kelvin(0.0)).0
        );
        assert_eq!(
            constants::SEA_LEVEL_TEMPERATURE.0 - 13.0,
            calculate_isa_temperature(Metres(2000.0), Kelvin(0.0)).0
        );
        assert_eq!(
            constants::TROPOPAUSE_TEMPERATURE.0,
            calculate_isa_temperature(constants::TROPOPAUSE_ALTITUDE, Kelvin(0.0)).0
        );
        assert!(
            libm::fabs(
                constants::TROPOPAUSE_TEMPERATURE.0 + 10.
                    - calculate_isa_temperature(constants::TROPOPAUSE_ALTITUDE, Kelvin(10.0)).0
            ) < 1.0e-9
        );
        assert_eq!(
            constants::TROPOPAUSE_TEMPERATURE.0,
            calculate_isa_temperature(Metres(12000.0), Kelvin(-10.0)).0
        );
    }

    #[test]
    fn test_calculate_density() {
        assert!(
            libm::fabs(
                constants::SEA_LEVEL_DENSITY.0
                    - calculate_density(
                        constants::SEA_LEVEL_PRESSURE,
                        constants::SEA_LEVEL_TEMPERATURE
                    )
                    .0
            ) < 2.0e-8
        );
        assert!(
            libm::fabs(
                0.3639176
                    - calculate_density(TROPOPAUSE_PRESSURE, constants::TROPOPAUSE_TEMPERATURE).0
            ) < 1.0e-6
        );
    }

    #[test]
    fn test_calculate_true_air_speed() {
        assert!(
            libm::fabs(
                150.0
                    - calculate_true_air_speed(
                        MetresPerSecond(150.0),
                        constants::SEA_LEVEL_PRESSURE,
                        constants::SEA_LEVEL_TEMPERATURE
                    )
                    .0
            ) < 1.0e-9
        );
        assert!(
            libm::fabs(
                164.458
                    - calculate_true_air_speed(
                        MetresPerSecond(150.0),
                        Pascals(79495.201),
                        Kelvin(constants::SEA_LEVEL_TEMPERATURE.0 - 13.0)
                    )
                    .0
            ) < 0.001
        );
    }

    #[test]
    fn test_calculate_calibrated_air_speed() {
        assert!(
            libm::fabs(
                150.0
                    - calculate_calibrated_air_speed(
                        MetresPerSecond(150.0),
                        constants::SEA_LEVEL_PRESSURE,
                        constants::SEA_LEVEL_TEMPERATURE
                    )
                    .0
            ) < 1.0e-9
        );
        assert!(
            libm::fabs(
                150.0
                    - calculate_calibrated_air_speed(
                        MetresPerSecond(164.458),
                        Pascals(79495.201),
                        Kelvin(constants::SEA_LEVEL_TEMPERATURE.0 - 13.0)
                    )
                    .0
            ) < 0.001
        );
    }

    #[test]
    fn test_speed_of_sound() {
        assert_eq!(0.0, speed_of_sound(Kelvin(0.0)).0);
        assert!(
            libm::fabs(
                constants::SEA_LEVEL_SPEED_OF_SOUND.0
                    - speed_of_sound(constants::SEA_LEVEL_TEMPERATURE).0
            ) < 0.001
        );
        assert!(libm::fabs(295.070 - speed_of_sound(constants::TROPOPAUSE_TEMPERATURE).0) < 0.001);
    }

    #[test]
    fn test_mach_true_air_speed() {
        assert!(
            libm::fabs(
                0.8 * constants::SEA_LEVEL_SPEED_OF_SOUND.0
                    - mach_true_air_speed(0.8, constants::SEA_LEVEL_TEMPERATURE).0
            ) < 0.001
        );
        assert!(
            libm::fabs(250.809 - mach_true_air_speed(0.85, constants::TROPOPAUSE_TEMPERATURE).0)
                < 0.001
        );
    }

    #[test]
    fn test_calculate_crossover_altitude() {
        let cas = MetresPerSecond(155.0);
        let crossover_altitude = calculate_crossover_altitude(cas, 0.79);
        assert!(libm::fabs(9070.814 - crossover_altitude.0) < 0.001);

        // The TAS should be the same from both CAS and MACH at the crossover_altitude
        let pressure = calculate_isa_pressure(crossover_altitude);
        let temperature = calculate_isa_temperature(crossover_altitude, Kelvin(0.0));
        let tas_from_cas = calculate_true_air_speed(cas, pressure, temperature);
        let tas_from_mach = mach_true_air_speed(0.79, temperature);
        assert!(libm::fabs(tas_from_cas.0 - tas_from_mach.0) < 0.001);
    }
}
