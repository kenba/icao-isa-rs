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

//! The ICAO Standard Atmosphere primary constants and characteristics.  
//! See Manual of the ICAO Standard Atmosphere; ICAO Doc 7488/3.

use icao_units::si::{
    Kelvin, KilogramsPerCubicMetre, Metres, MetresPerSecond, MetresPerSecondSquared, Pascals,
};

// Constants from ICAO Doc 7488/3, Table A.

/// The acceleration due to gravity at latitude 45°32'33'' using Lambert's equation.
#[allow(non_upper_case_globals)]
pub const g: MetresPerSecondSquared = MetresPerSecondSquared(9.806_65);

/// The adiabatic index of air (K), dimensionless.
pub const K: f64 = 1.4;

/// The real gas constant for air (R) in metres squared per Kelvin seconds squared.
pub const R: f64 = 287.052_87;

/// ISA Sea level temperature.
pub const SEA_LEVEL_TEMPERATURE: Kelvin = Kelvin(288.15);

/// ISA Sea level pressure.
pub const SEA_LEVEL_PRESSURE: Pascals = Pascals(101_325.0);

/// ISA Sea level density.
pub const SEA_LEVEL_DENSITY: KilogramsPerCubicMetre = KilogramsPerCubicMetre(1.225);

/// ISA Sea level speed of sound (a0).
/// ICAO Doc 7488/3, Table C.
pub const SEA_LEVEL_SPEED_OF_SOUND: MetresPerSecond = MetresPerSecond(340.294);

// Constants from ICAO Doc 7488/3, Table D.

/// ISA tropopause temperature.
pub const TROPOPAUSE_TEMPERATURE: Kelvin = Kelvin(216.65);

/// The ISA temperature gradient from Sea level to the tropopause altitude in K/m.  
/// AKA Lapse Rate.
pub const TEMPERATURE_GRADIENT: f64 = -0.0065;

/// ISA tropopause altitude.
pub const TROPOPAUSE_ALTITUDE: Metres = Metres(11000.0);
