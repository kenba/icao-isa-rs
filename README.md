# icao-isa

[![crates.io](https://img.shields.io/crates/v/icao-isa.svg)](https://crates.io/crates/icao-isa)
[![docs.io](https://docs.rs/icao-isa/badge.svg)](https://docs.rs/icao-isa/)
[![License](https://img.shields.io/badge/License-MIT-blue)](https://opensource.org/license/mit/)
[![Rust](https://github.com/kenba/icao-isa-rs/actions/workflows/rust.yml/badge.svg)](https://github.com/kenba/icao-isa-rs/actions)
[![codecov](https://codecov.io/gh/kenba/icao-isa-rs/graph/badge.svg?token=6DTOY9Y4BT)](https://codecov.io/gh/kenba/icao-isa-rs)

An implementation of the [International Civil Aviation Organization](https://icao.int/) (ICAO) [International Standard Atmosphere](https://en.wikipedia.org/wiki/International_Standard_Atmosphere) (ISA), see [ICAO Doc 7488/3](https://standart.aero/en/icao/book/doc-7488-manual-of-the-icao-standard-atmosphere-extended-to-80-kilometres-262-500-feet-en-cons).

The library also includes functions for calculating:

- true airspeed ([TAS](https://en.wikipedia.org/wiki/True_airspeed)) from calibrated airspeed ([CAS](https://en.wikipedia.org/wiki/Calibrated_airspeed)), pressure and temperature;
- CAS from TAS, pressure and temperature;
- TAS from [Mach number](https://en.wikipedia.org/wiki/Mach_number) and temperature;
- and the crossover altitude between CAS / Mach flight regimes.

The equations for the functions above are from [BADA User Manual revision 3-12](https://www.scribd.com/document/289480324/1-User-Manual-Bada-3-12).

The library is declared [no_std](https://docs.rust-embedded.org/book/intro/no-std.html)
so it can be used in embedded applications.

## Contribution

If you want to contribute through code or documentation, the [Contributing](CONTRIBUTING.md) guide is the best place to start.  
Just please abide by our [Code of Conduct](CODE_OF_CONDUCT.md).

## License

`icao-isa-rs` is provided under a MIT license, see [LICENSE](LICENSE).
