use std::time::{SystemTime, UNIX_EPOCH};

fn main() {
    // Get the timestamp in seconds
    let now = SystemTime::now();
    let duration_since_epoch = now.duration_since(UNIX_EPOCH).expect("Time went backwards");
    let timestamp_in_seconds = duration_since_epoch.as_secs();

    // get the list of phases surrounding a date 
    let phase_max_list = moon_calc::phasehunt(timestamp_in_seconds as f64);
    println!("{:?}", phase_max_list);
}
