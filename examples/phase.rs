use std::time::{SystemTime, UNIX_EPOCH};

fn main() {
    // Get the current time
    let now = SystemTime::now();
    let duration_since_epoch = now.duration_since(UNIX_EPOCH).expect("Time went backwards");
    let timestamp_in_seconds = duration_since_epoch.as_secs();

    println!("Current UTC timestamp in seconds: {}", timestamp_in_seconds);

    let x = moon_calc::phase(timestamp_in_seconds as f64);
    println!("{:?}", x.get_icon());
    println!("{:?}", x.get_name());
    println!("{:?}", x);
}
