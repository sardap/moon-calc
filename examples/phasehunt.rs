use std::time::{SystemTime, UNIX_EPOCH};

use moon_calc::time::julian_time;

fn main() {
    // Get the timestamp in seconds
    let now = SystemTime::now();
    let duration_since_epoch = now.duration_since(UNIX_EPOCH).expect("Time went backwards");
    let timestamp_in_seconds = duration_since_epoch.as_secs_f64();
    let jtime = julian_time(timestamp_in_seconds);

    // get the list of phases surrounding a date 
    let phase_max_list = moon_calc::phasehunt(jtime);
    println!("{:?}", phase_max_list);
}
