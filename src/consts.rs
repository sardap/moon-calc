use crate::Degree;
pub use std::f64::consts::PI;

/// 1980 January 0.0
pub const EPOCH: f64 = 2444238.5;

/// 1970 January 1.1
pub const JULIAN_DATE_EPOCH: f64 = 2440587.5;

/*
    Constants defining the Sun's apparent orbit.
*/
/// Ecliptic longitude of the Sun at epoch 1980.0
pub const SUN_ECLIPTIC_LONGITUDE_AT_EPOCH: f64 = 278.833540;
/// Ecliptic longitude of the Sun at perigee
pub const SUN_ECLIPTIC_LONGITUDE_AT_PERIGEE: f64 = 282.596403;
/// Eccentricity of Earth's orbit
pub const ECCENTRICITY: f64 = 0.016718;
/// Semi-major axis of Earth's orbit, km
pub const SUN_SEMI_MAJOR_AXIS: f64 = 1.495985e8;
/// Sun's angular size, degrees, at semi-major axis distance
pub const SUN_ANGULAR_SIZE: Degree = 0.533128;

/*
    Elements of the Moon's orbit, epoch 1980.0.
*/
/// moon's mean longitude at the epoch
pub const MOON_ECLIPTIC_LONGITUDE_AT_EPOCH: f64 = 64.975464;
/// mean longitude of the perigee at the epoch
pub const MOON_ECLIPTIC_LONGITUDE_AT_PERIGEE: f64 = 349.383063;
/// mean longitude of the node at the epoch
pub const MEAN_LONGITUDE_NODE: f64 = 151.950429;
/// inclination of the Moon's orbit
pub const MOON_INCLINATION: f64 = 5.145396;
/// eccentricity of the Moon's orbit
pub const MOON_ECCENTRICITY: f64 = 0.054900;
/// moon's angular size at distance a from Earth
pub const MOON_ANGULAR_SIZE: f64 = 0.5181;
/// semi-major axis of Moon's orbit in km
pub const MOON_SEMI_MAJOR_AXIS: f64 = 384401.0;
/// parallax at distance a from Earth
pub const MOON_PARALLAX: f64 = 0.9507;
/// synodic month (new Moon to new Moon) in days
pub const SYNODIC_MONTH: f64 = 29.53058868;

/// Name of the phase
pub const PHASE_NAME: [&str; 8] = [
    "New Moon",        // 0
    "Waxing Crescent", // 1
    "First Quarter",   // 2
    "Waxing Gibbous",  // 3
    "Full Moon",       // 4
    "Waning Gibbous",  // 3
    "Last Quarter",    // 2
    "Waning Crescent", // 1
];

/// Unicode representation of the given phase
pub const MOON_ICON: [&str; 8] = [
    "\u{1f311}",
    "\u{1f312}",
    "\u{1f313}",
    "\u{1f314}",
    "\u{1f315}",
    "\u{1f316}",
    "\u{1f317}",
    "\u{1f318}",
];
