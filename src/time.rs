use crate::JULIAN_DATE_EPOCH;

/// converts the given UNIX timestamp into the corresponding Julian date.
///
/// # Parameters:
/// - `timestamp`: The time in seconds since the UNIX epoch (January 1, 1970, 00:00:00 UTC) as a `f64` value.
///
/// # Returns:
/// - A `f64` representing the Julian date with day fraction, which is the number of days since January 1, 4713 BCE, 12:00 UTC,
///   as used in astronomical calculations.
pub fn julian_time(timestamp: f64) -> f64 {
    (timestamp / 86400.0_f64) + JULIAN_DATE_EPOCH
}

/// converts a Julian date (including the fractional day) to the equivalent UNIX timestamp
///
/// # Parameters:
/// - `julian_date`: The Julian date to convert, as a `f64`. This includes both the date and the fractional part representing the time of day.
///
/// # Returns:
/// - A `f64` representing the number of seconds elapsed since the UNIX epoch (January 1, 1970, 00:00:00 UTC).
pub fn julian_days_to_unix_seconds(julian_date: f64) -> f64 {
    (julian_date - JULIAN_DATE_EPOCH) * 86400.0
}

/// convert Julian date to year, month, day
// for the seemingly random values read: https://aa.usno.navy.mil/faq/JD_formula
// POSSIBLE ISSUE the return type on perl is float in c its int see https://play.rust-lang.org/?version=stable&mode=debug&edition=2021&gist=35e90c45e56ab65b8c0d5bc10cee42d7
pub fn julian_year(jtimestamp: f64) -> (i32, i32, f64) {
    let jtimestamp: f64 = jtimestamp + 0.5;
    let int_part: f64 = f64::floor(jtimestamp);
    let float_part: f64 = jtimestamp - int_part;

    // if the date is before the Gregorian calendar transition date (October 15, 1582)
    let a = if int_part < 2299161.0_f64 {
        int_part
    } else {
        // Calculate the correction factor for the Gregorian calendar
        let alpha = f64::floor((int_part - 1867216.25) / 36524.25);
        int_part + 1.0 + alpha - f64::floor(alpha / 4.0)
    };

    let b = a + 1524.0;
    let c = f64::floor((b - 122.1) / 365.25);
    let d = f64::floor(365.25 * c);
    let e = f64::floor((b - d) / 30.6001);

    let day = b - d - f64::floor(30.6001 * e) + float_part;
    let month = if e < 14.0 { e - 1.0 } else { e - 13.0 };
    let year = if month > 2.0 { c - 4716.0 } else { c - 4715.0 };

    (year as i32, month as i32, day)
}
