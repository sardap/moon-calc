use std::time::{Duration, SystemTime, UNIX_EPOCH};

use moon_calc::time::julian_time;

fn main() {
    // Get the current time
    let now: SystemTime = SystemTime::now();
    let now_sec = now.duration_since(UNIX_EPOCH).unwrap().as_secs_f64();
    let now_jtime = julian_time(now_sec);

    let date2 = UNIX_EPOCH + Duration::from_secs(100);
    let date2_sec = date2.duration_since(UNIX_EPOCH).unwrap().as_secs_f64();
    let date2_jtime = julian_time(date2_sec);


    let phaselist = moon_calc::phaselist(date2_jtime, now_jtime);

    println!("{:?}", phaselist);
}
