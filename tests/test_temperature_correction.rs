// Copyright (c) 2026 Via Technology Ltd.

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

//! Altitude temoperature correction tests using values from:
//! [Eurocae ED-323 Minimum Operational Performance Standards - Required Navigation Performance for Area Navigation](https://www.eurocae.net/product/ed-323-minimum-operational-performance-standards-for-required-navigation-performance-for-area-navigation/)
//! Appendix H Table H-2: Altitude Correction Values

extern crate icao_isa;

use icao_isa::{Kelvin, Metres, calculate_temperature_correction_delta_altitude};
use icao_units::non_si::Feet;

const TEST_ALTITUDES_FT: [f64; 5] = [250.0, 1500.0, 2500.0, 5000.0, 10000.0];
const CORRECTIONS_P25_0_FT: [i32; 5] = [-20, -120, -201, -405, -822];
const CORRECTIONS_P10_0_FT: [i32; 5] = [-8, -51, -85, -170, -347];
const CORRECTIONS_M15_0_FT: [i32; 5] = [14, 83, 139, 280, 571];
const CORRECTIONS_M30_0_FT: [i32; 5] = [29, 175, 294, 594, 1215];
const CORRECTIONS_M45_0_FT: [i32; 5] = [46, 279 + 1, 468, 948, 1947]; // 1500ft value 1 higher
const CORRECTIONS_M60_0_FT: [i32; 5] = [66, 398, 667, 1352, 2787];

const CORRECTIONS_P35_5000_FT: [i32; 5] = [-28, -168, -281, -565, -1147];
const CORRECTIONS_P20_5000_FT: [i32; 5] = [-17, -101, -168, -339, -689];
const CORRECTIONS_P10_5000_FT: [i32; 5] = [-9, -52, -87, -175, -355 - 1]; // 10000ft value 1 lower
const CORRECTIONS_M5_5000_FT: [i32; 5] = [5, 28, 47, 95, 194];
const CORRECTIONS_M20_5000_FT: [i32; 5] = [19, 118, 197, 398, 813];
const CORRECTIONS_M35_5000_FT: [i32; 5] = [36, 218, 365, 739, 1516];
const CORRECTIONS_M50_5000_FT: [i32; 5] = [55, 332, 556, 1128, 2322];

const CORRECTIONS_P55_15000_FT: [i32; 5] = [-44, -263, -440, -885, -1794];
const CORRECTIONS_P40_15000_FT: [i32; 5] = [-33, -201, -335, -676, -1371];
const CORRECTIONS_P30_15000_FT: [i32; 5] = [-26, -155, -260, -524, -1064];
const CORRECTIONS_P15_15000_FT: [i32; 5] = [-13, -81, -136, -274, -557 - 1]; // 10000ft value 1 lower
const CORRECTIONS_M03_15000_FT: [i32; 5] = [0, 2, 3, 6, 12];
const CORRECTIONS_M15_15000_FT: [i32; 5] = [16, 95, 159, 322, 658];
const CORRECTIONS_M30_15000_FT: [i32; 5] = [33, 201, 336, 681, 1398];
const CORRECTIONS_M45_15000_FT: [i32; 5] = [53, 321, 539, 1094, 2256];

fn check_temperature_corrections(
    delta_temperature: Kelvin<f64>,
    ref_elevation: Metres<f64>,
    corrections: &[i32],
) {
    let tolerance = 0.5;

    for (altitude, expected) in TEST_ALTITUDES_FT.into_iter().zip(corrections) {
        let result = calculate_temperature_correction_delta_altitude(
            Metres::from(Feet(altitude)) + ref_elevation,
            delta_temperature,
            ref_elevation,
            Metres::from(Feet(tolerance)),
        );
        let expected = *expected as f64;
        let error = (Feet::from(result.0).0 - expected).abs();
        assert!(error < tolerance)
    }
}

#[test]
fn test_temperature_correction_sea_level() {
    // Reference Sea Level
    check_temperature_corrections(Kelvin(25.0), Metres(0.0), &CORRECTIONS_P25_0_FT);
    check_temperature_corrections(Kelvin(10.0), Metres(0.0), &CORRECTIONS_P10_0_FT);
    check_temperature_corrections(Kelvin(-15.0), Metres(0.0), &CORRECTIONS_M15_0_FT);
    check_temperature_corrections(Kelvin(-30.0), Metres(0.0), &CORRECTIONS_M30_0_FT);
    check_temperature_corrections(Kelvin(-45.0), Metres(0.0), &CORRECTIONS_M45_0_FT);
    check_temperature_corrections(Kelvin(-60.0), Metres(0.0), &CORRECTIONS_M60_0_FT);
}

#[test]
fn test_temperature_correction_from_5000ft() {
    // Reference 5000 ft
    let feet_5000 = Metres::from(Feet(5000.0));
    check_temperature_corrections(Kelvin(34.9), feet_5000, &CORRECTIONS_P35_5000_FT);
    check_temperature_corrections(Kelvin(19.9), feet_5000, &CORRECTIONS_P20_5000_FT);
    check_temperature_corrections(Kelvin(9.9), feet_5000, &CORRECTIONS_P10_5000_FT);
    check_temperature_corrections(Kelvin(-5.1), feet_5000, &CORRECTIONS_M5_5000_FT);
    check_temperature_corrections(Kelvin(-20.1), feet_5000, &CORRECTIONS_M20_5000_FT);
    check_temperature_corrections(Kelvin(-35.1), feet_5000, &CORRECTIONS_M35_5000_FT);
    check_temperature_corrections(Kelvin(-50.1), feet_5000, &CORRECTIONS_M50_5000_FT);
}

#[test]
fn test_temperature_correction_from_15000ft() {
    // Reference 15000 ft
    let feet_15000 = Metres::from(Feet(15000.0));
    check_temperature_corrections(Kelvin(54.7), feet_15000, &CORRECTIONS_P55_15000_FT);
    check_temperature_corrections(Kelvin(39.7), feet_15000, &CORRECTIONS_P40_15000_FT);
    check_temperature_corrections(Kelvin(29.7), feet_15000, &CORRECTIONS_P30_15000_FT);
    check_temperature_corrections(Kelvin(14.7), feet_15000, &CORRECTIONS_P15_15000_FT);
    check_temperature_corrections(Kelvin(-0.3), feet_15000, &CORRECTIONS_M03_15000_FT);
    check_temperature_corrections(Kelvin(-15.3), feet_15000, &CORRECTIONS_M15_15000_FT);
    check_temperature_corrections(Kelvin(-30.3), feet_15000, &CORRECTIONS_M30_15000_FT);
    check_temperature_corrections(Kelvin(-45.3), feet_15000, &CORRECTIONS_M45_15000_FT);
}
