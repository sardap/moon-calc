use std::time::{Duration, SystemTime, UNIX_EPOCH};

fn main() {
    // Get the current time
    let now: SystemTime = SystemTime::now();
    let now_sec = now.duration_since(UNIX_EPOCH).unwrap().as_secs_f64();

    let date2 = UNIX_EPOCH + Duration::from_secs(100);
    let date2_sec = date2.duration_since(UNIX_EPOCH).unwrap().as_secs_f64();


    let phaselist = moon_calc::phaselist(date2_sec, now_sec);

    println!("{:?}", phaselist);
}
