#[warn(missing_docs)]
/// constants
pub mod consts;

/// utilities for time
pub mod time;

use consts::*;
use time::*;

pub type Degree = f64;
pub type Radiant = f64;

/// deg to rad
fn to_rad(val: Degree) -> Radiant {
    val * (PI / 180.0)
}
/// rad to deg
fn to_deg(val: Radiant) -> Degree {
    val * (180.0 / PI)
}

/// sin from degree
fn dsin(val: Degree) -> f64 {
    to_rad(val).sin()
}
/// cos from degree
fn dcos(val: Degree) -> f64 {
    to_rad(val).cos()
}

/// fix angle when it lapses
fn fix_angle(angle: f64) -> f64 {
    angle - 360.0 * (angle / 360.0).floor()
}

/// Represents information about the moon, including its julian date, phase,
/// age, illumination, distance, and lunation.
#[derive(Debug, Clone, Copy)]
pub struct Moon {
    /// A continuous count of days and fractions since noon Universal Time on January 1, 4713 BC
    pub julian_date: f64,
    /// Phase of the moon.
    pub phase: f64,
    /// Age of the moon.
    pub age: f64,
    /// Illumination of the moon (0 to 1 where 1 is a full moon).
    pub illumination: f64,
    /// Distance of the moon in earth radii.
    pub distance: f64,
    /// angle of the moon
    pub angle: f64,
    /// lunation number.
    pub lunation: u16,
}

/// Kepler's equation for eccentric anomaly (e) given the mean anomaly (M)
/// and the orbital eccentricity.
pub fn kepler(mean_anomaly: f64, eccentricity: f64) -> f64 {
    let mean_anomaly = to_rad(mean_anomaly);
    let mut eccentric_anomaly = mean_anomaly; // as a guess
    const EPSILON: f64 = 1e-6; // tolerance

    // solving Kepler's equation,
    // using: https://en.wikipedia.org/wiki/Kepler%27s_equation#Newton%27s_method
    while {
        let delta = eccentric_anomaly - eccentricity * eccentric_anomaly.sin() - mean_anomaly;
        eccentric_anomaly -= delta / (1.0 - eccentricity * eccentric_anomaly.cos());

        delta.abs() > EPSILON
    } {}
    eccentric_anomaly
}

/// Calculate phase of moon as a fraction:
/// # Parameters
/// * `date_timestamp` - Time of the state returned, in julian date and fraction
///
/// Returns the terminator phase angle as a percentage of a full circle (i.e., 0 to 1),
/// and stores into pointer arguments the illuminated fraction of
/// the Moon's disc, the Moon's age in days and fraction, the
/// distance of the Moon from the centre of the Earth, and the
/// angular diameter subtended by the Moon as seen by an observer
/// at the centre of the Earth.
pub fn phase(date_timestamp: f64) -> PhaseRes {
    // date within epoch
    let jd = date_timestamp - EPOCH;

    // Calculate the Sun's position.
    let sun_mean_anomaly = fix_angle((360.0 / 365.2422) * jd); // mean anomaly of the Sun

    // converted from perigee
    // co-ordinates to epoch 1980.0
    let adjusted_mean_anomaly = fix_angle(
        sun_mean_anomaly + SUN_ECLIPTIC_LONGITUDE_AT_EPOCH - SUN_ECLIPTIC_LONGITUDE_AT_PERIGEE,
    );

    // true anomaly
    let true_anomaly: Degree = {
        let ec = kepler(adjusted_mean_anomaly, ECCENTRICITY); // solve equation of Kepler
        let ec = ((1.0 + ECCENTRICITY) / (1.0 - ECCENTRICITY)).sqrt() * (ec / 2.0).tan();

        2.0 * to_deg(ec.atan())
    };

    // Sun's geocentric ecliptic longitude
    let lambda_sun = fix_angle(true_anomaly + SUN_ECLIPTIC_LONGITUDE_AT_PERIGEE);

    // Orbital distance factor
    let orbital_distance_factor =
        (1.0 + ECCENTRICITY * dcos(true_anomaly)) / (1.0 - ECCENTRICITY * ECCENTRICITY);
    let sun_dist = SUN_SEMI_MAJOR_AXIS / orbital_distance_factor; // distance to Sun in km
    let sun_ang: Degree = orbital_distance_factor * SUN_ANGULAR_SIZE; // Sun's angular size in degrees

    // Calculation of the Moon's position

    // Moon's mean longitude
    let moon_mean_longitude = fix_angle(13.1763966 * jd + MOON_ECLIPTIC_LONGITUDE_AT_EPOCH);
    // Moon's mean anomaly
    let moon_mean_anomaly =
        fix_angle(moon_mean_longitude - 0.1114041 * jd - MOON_ECLIPTIC_LONGITUDE_AT_PERIGEE);

    // Evection
    let moon_evection = 1.2739 * dsin(2.0 * (moon_mean_longitude - lambda_sun) - moon_mean_anomaly);

    // Annual equation
    let annual_equation = 0.1858 * dsin(adjusted_mean_anomaly);

    // Correction term
    let a3 = 0.37 * dsin(adjusted_mean_anomaly);

    // Corrected anomaly
    let mmp = moon_mean_anomaly + moon_evection - annual_equation - a3;

    // Correction for the equation of the centre
    let mec = 6.2886 * dsin(mmp);

    // Another correction term
    let a4 = 0.214 * dsin(2.0 * mmp);

    // True longitude
    let true_longitude = {
        // Corrected longitude
        let corrected_longitude = moon_mean_longitude + moon_evection + mec - annual_equation + a4;
        // Variation.
        let variation = 0.6583 * dsin(2.0 * (corrected_longitude - lambda_sun));

        corrected_longitude + variation
    };

    // Moon's ascending node mean longitude
    let node_mean_longitude = fix_angle(MEAN_LONGITUDE_NODE - 0.0529539 * jd);

    // Corrected longitude of the node
    let node_corrected_longitude = node_mean_longitude - 0.16 * dsin(adjusted_mean_anomaly);

    let _lambda_moon: Degree = {
        // Y inclination coordinate
        let y: f64 = dsin(true_longitude - node_corrected_longitude) * dcos(MOON_INCLINATION);
        // X inclination coordinate
        let x: f64 = dcos(true_longitude - node_corrected_longitude);

        // Ecliptic longitude.
        to_deg(f64::atan2(y, x) + node_corrected_longitude)
    };

    // Ecliptic latitude
    let _beta_m = to_deg(f64::asin(
        dsin(true_longitude - node_corrected_longitude) * dsin(MOON_INCLINATION),
    ));

    // Calculation of the phase of the Moon

    // Age of the Moon in degrees
    let moon_age: Degree = true_longitude - lambda_sun;

    // Phase of the Moon
    let moon_phase = (1.0 - dcos(moon_age)) / 2.0;

    // Calculate distance of moon from the centre of the Earth
    let moon_dist = (MOON_SEMI_MAJOR_AXIS * (1.0 - MOON_ECCENTRICITY * MOON_ECCENTRICITY))
        / (1.0 + MOON_ECCENTRICITY * dcos(mmp + mec));

    // Calculate Moon's angular diameter
    let moon_diameter_fract = moon_dist / MOON_SEMI_MAJOR_AXIS;
    let angular_diameter = MOON_ANGULAR_SIZE / moon_diameter_fract;

    // Calculate Moon's parallax
    let _moon_parallax = MOON_PARALLAX / moon_diameter_fract;

    PhaseRes {
        phase: moon_phase,
        illuminated_fraction: fix_angle(moon_age) / 360.0,
        age: SYNODIC_MONTH * (fix_angle(moon_age) / 360.0),
        distance: moon_dist,
        angular_diameter,
        sun_distance: sun_dist,
        sun_angular_diameter: sun_ang,
    }
}

#[derive(Debug)]
#[non_exhaustive]
pub struct PhaseRes {
    /// Terminator phase angle as a percentage of a full circle (i.e., 0 to 1)
    pub phase: f64,
    /// Illuminated fraction of the Moon's disc
    pub illuminated_fraction: f64,
    /// Age of moon in days and fractions
    pub age: f64,
    /// Distance of the Moon from the center of the Earth in km
    pub distance: f64,
    /// Angular diameter in degrees
    pub angular_diameter: f64,
    /// Distance to Sun in km
    pub sun_distance: f64,
    /// Sun's angular diameter
    pub sun_angular_diameter: f64,
}

/// Calculates the time of the mean new Moon for a given base date.
///
/// The function computes the mean new Moon time based on a precomputed
/// synodic month index `k`, where `k = (year - 1900) * 12.3685`, and the
/// base date `sdate` in julian datw format.
///
/// # Parameters
/// - `sdate`: The base date in Julian date format.
/// - `k`: The precomputed synodic month index.
///
/// # Returns
/// The time of the mean new Moon in julian date format.
pub fn meanphase(sdate: f64, k: f64) -> f64 {
    // Time in julian centuries from 1900 January 0.5
    let t: f64 = (sdate - 2415020.0) / 36525.0;
    let t2 = t * t;
    let t3 = t2 * t;

    2415020.75933 + SYNODIC_MONTH * k + 0.0001178 * t2 - 0.000000155 * t3
        + 0.00033 * dsin(166.56 + 132.87 * t - 0.009173 * t2)
}

/// Enum representing the different phase selectors for the Moon.
///
/// The phase values correspond to:
/// - 0.0:  New Moon
/// - 0.25: First Quarter
/// - 0.5:  Full Moon
/// - 0.75: Last Quarter
#[derive(Debug, Clone, Copy)]
pub enum PhaseSelector {
    NewMoon,
    FirstQuarter,
    FullMoon,
    LastQuarter,
}
impl From<f64> for PhaseSelector {
    fn from(value: f64) -> Self {
        match value {
            0.0..=0.25 => Self::NewMoon,
            0.25..=0.50 => Self::FirstQuarter,
            0.50..=0.75 => Self::FullMoon,
            0.75..=1.0 => Self::LastQuarter,
            _ => panic!("phase cant be over 1.0 or under 0.0"),
        }
    }
}

impl PhaseSelector {
    /// Converts a `PhaseSelector` enum variant to its corresponding phase value (f64).
    fn to_value(self) -> f64 {
        match self {
            PhaseSelector::NewMoon => 0.0,
            PhaseSelector::FirstQuarter => 0.25,
            PhaseSelector::FullMoon => 0.5,
            PhaseSelector::LastQuarter => 0.75,
        }
    }

    /// A vec of all possible values
    #[allow(dead_code)]
    fn all() -> Vec<Self> {
        vec![
            PhaseSelector::NewMoon,
            PhaseSelector::FirstQuarter,
            PhaseSelector::FullMoon,
            PhaseSelector::LastQuarter,
        ]
    }
}

/// Calculates the true phase time of the Moon based on the mean new Moon time and phase selector.
///
/// Given a `k` value (used to determine the mean phase of the new Moon) and a phase selector
/// (0.0 for new moon, 0.25 for first quarter, 0.5 for full moon, 0.75 for last quarter),
/// this function computes the corrected phase time.
///
/// # Parameters
/// - `k`: The precomputed synodic month index with phase added.
/// - `phase`: A phase selector (0.0, 0.25, 0.5, or 0.75).
///
/// # Returns
/// The true corrected phase time in Julian date format.
pub fn truephase(mut k: f64, phase: PhaseSelector) -> f64 {
    let phase = phase.to_value();

    k += phase; // Add phase to new moon time
    let t = k / 1236.85; // time in Julian centuries from 1900 January 0.5

    let (mut mean_phase_time, sun_mean_anomaly, moon_mean_anomaly, moon_argument_latitude) = {
        let t2 = t * t;
        let t3 = t2 * t;

        // Mean time of phases
        let mean_phase_time = 2415020.75933 + SYNODIC_MONTH * k + 0.0001178 * t2 - 0.000000155 * t3
            + 0.00033 * dsin(166.56 + 132.87 * t - 0.009173 * t2);

        // Sun's mean anomaly
        let sun_mean_anomaly = 359.2242 + 29.10535608 * k - 0.0000333 * t2 - 0.00000347 * t3;

        // Moon's mean anomaly
        let moon_mean_anomaly = 306.0253 + 385.81691806 * k + 0.0107306 * t2 + 0.00001236 * t3;

        // Moon's argument of latitude
        let moon_argument_latitude = 21.2964 + 390.67050646 * k - 0.0016528 * t2 - 0.00000239 * t3;

        (
            mean_phase_time,
            sun_mean_anomaly,
            moon_mean_anomaly,
            moon_argument_latitude,
        )
    };
    // Correcctions
    if (phase < 0.01) || ((phase - 0.5).abs() < 0.01) {
        // Corrections for New and Full Moon.
        mean_phase_time += (0.1734 - 0.000393 * t) * dsin(sun_mean_anomaly)
            + 0.0021 * dsin(2.0 * sun_mean_anomaly)
            - 0.4068 * dsin(moon_mean_anomaly)
            + 0.0161 * dsin(2.0 * moon_mean_anomaly)
            - 0.0004 * dsin(3.0 * moon_mean_anomaly)
            + 0.0104 * dsin(2.0 * moon_argument_latitude)
            - 0.0051 * dsin(sun_mean_anomaly + moon_mean_anomaly)
            - 0.0074 * dsin(sun_mean_anomaly - moon_mean_anomaly)
            + 0.0004 * dsin(2.0 * moon_argument_latitude + sun_mean_anomaly)
            - 0.0004 * dsin(2.0 * moon_argument_latitude - sun_mean_anomaly)
            - 0.0006 * dsin(2.0 * moon_argument_latitude + moon_mean_anomaly)
            + 0.0010 * dsin(2.0 * moon_argument_latitude - moon_mean_anomaly)
            + 0.0005 * dsin(sun_mean_anomaly + 2.0 * moon_mean_anomaly);
    } else if (phase - 0.25).abs() < 0.01 || ((phase - 0.75).abs() < 0.01) {
        mean_phase_time += (0.1721 - 0.0004 * t) * dsin(sun_mean_anomaly)
            + 0.0021 * dsin(2.0 * sun_mean_anomaly)
            - 0.6280 * dsin(moon_mean_anomaly)
            + 0.0089 * dsin(2.0 * moon_mean_anomaly)
            - 0.0004 * dsin(3.0 * moon_mean_anomaly)
            + 0.0079 * dsin(2.0 * moon_argument_latitude)
            - 0.0119 * dsin(sun_mean_anomaly + moon_mean_anomaly)
            - 0.0047 * dsin(sun_mean_anomaly - moon_mean_anomaly)
            + 0.0003 * dsin(2.0 * moon_argument_latitude + sun_mean_anomaly)
            - 0.0004 * dsin(2.0 * moon_argument_latitude - sun_mean_anomaly)
            - 0.0006 * dsin(2.0 * moon_argument_latitude + moon_mean_anomaly)
            + 0.0021 * dsin(2.0 * moon_argument_latitude - moon_mean_anomaly)
            + 0.0003 * dsin(sun_mean_anomaly + 2.0 * moon_mean_anomaly)
            + 0.0004 * dsin(sun_mean_anomaly - 2.0 * moon_mean_anomaly)
            - 0.0003 * dsin(2.0 * sun_mean_anomaly + moon_mean_anomaly);
        if phase < 0.5 {
            //  First quarter correction.
            mean_phase_time +=
                0.0028 - 0.0004 * dcos(sun_mean_anomaly) + 0.0003 * dcos(moon_mean_anomaly);
        } else {
            // Last quarter correction.
            mean_phase_time +=
                -0.0028 + 0.0004 * dcos(sun_mean_anomaly) - 0.0003 * dcos(moon_mean_anomaly);
        }
    }

    mean_phase_time
}

/// phasehunt finds the time of phases of the moon which surround the current
/// date.
/// Five phases are found, starting and ending with the
/// new moons which bound the current lunation
pub fn phasehunt(starting_date: f64) -> (f64, f64, f64, f64, f64) {
    let adate = starting_date - 45.0;
    let (year, month, _day) = julian_year(adate);

    let mut k1 =
        f64::floor((year as f64 + ((month as f64 - 1.0) * (1.0 / 12.0)) - 1900.0) * 12.3685);
    let mut nt1 = meanphase(adate, k1);
    let (mut nt2, mut k2);

    let mut current_phase_date = nt1;
    loop {
        current_phase_date += SYNODIC_MONTH;
        k2 = k1 + 1.0;
        nt2 = meanphase(current_phase_date, k2);
        if (nt1 <= starting_date) && (nt2 > starting_date) {
            break;
        }
        nt1 = nt2;
        k1 = k2;
    }

    // this needs improvements
    (
        julian_days_to_unix_seconds(truephase(k1, PhaseSelector::NewMoon)),
        julian_days_to_unix_seconds(truephase(k1, PhaseSelector::FirstQuarter)),
        julian_days_to_unix_seconds(truephase(k1, PhaseSelector::FullMoon)),
        julian_days_to_unix_seconds(truephase(k1, PhaseSelector::LastQuarter)),
        julian_days_to_unix_seconds(truephase(k2, PhaseSelector::NewMoon)),
    )
}

/// phaselist finds the time of phases of the moon between two dates
/// # Parameters
/// * `start_date` -  the starting date of the range, in seconds since Unix epoch
/// * `end_date` - the ending date of the range, in seconds since Unix epoch
pub fn phaselist(start_date: f64, end_date: f64) -> Vec<f64> {
    let sd_sec = julian_time(start_date);
    let ed_sec = julian_time(end_date);

    let mut phases = Vec::new();
    let (year, month, _) = julian_year(sd_sec);

    let mut k =
        (((year as f64 + ((month as f64 - 1.0) * (1.0 / 12.0))) - 1900.0) * 12.3685) as i32 - 2;

    loop {
        k += 1;

        for &phase in PhaseSelector::all().iter() {
            let day: f64 = truephase(k as f64, phase);

            if day >= ed_sec {
                return phases;
            }

            if day >= sd_sec {
                if phases.is_empty() {
                    phases.push(4.0 * phase.to_value());
                }
                phases.push(julian_days_to_unix_seconds(day));
            }
        }
    }
}

impl PhaseRes {
    /// Get the icon of the phase
    pub fn get_icon(&self) -> &str {
        MOON_ICON
            .get(self.get_nth().unwrap())
            .expect("bad implementation the index is out of bounds.")
    }
    /// Get the name of the phase
    pub fn get_name(&self) -> &str {
        PHASE_NAME
            .get(self.get_nth().unwrap())
            .expect("bad implementation the index is out of bounds.")
    }

    // the fuk is dis
    // idk :p
    fn get_nth(&self) -> Result<usize, String> {
        let offset = (1. / 8.) / 2.;
        let phase_ranges = [
            ((0.0 - offset)..=(0.125 - offset), 0),
            ((0.125 - offset)..=(0.25 - offset), 1),
            ((0.25 - offset)..=(0.375 - offset), 2),
            ((0.375 - offset)..=(0.50 - offset), 3),
            ((0.50 - offset)..=(0.625 - offset), 4),
            ((0.625 - offset)..=(0.75 - offset), 5),
            ((0.75 - offset)..=(0.875 - offset), 6),
            ((0.875 - offset)..=(1.0 - offset), 7),
        ];

        for (r, n) in phase_ranges {
            if r.contains(&self.phase) {
                return Ok(n);
            }
        }
        Err("value must be between 0..=1".to_string())
    }
}
