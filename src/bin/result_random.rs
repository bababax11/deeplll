use std::path::Path;

use deeplll::parse::{cat, parse_sec_len_cnt, read_dir};
use once_cell::sync::Lazy;
use regex::Regex;

fn parse_path_to_tuple_random(path: &Path) -> (&str, u32, u64, u64, u32, &str) {
    let path_str: Vec<_> = path.to_str().unwrap().split('/').collect();
    let algo = path_str[1];
    let path_str = path_str[2];
    static RE_TUPLE: Lazy<Regex> = Lazy::new(|| Regex::new(r"[\.0-9]+").unwrap());

    let mut caps = RE_TUPLE.captures_iter(path_str);
    let (all_dim, seed, cnt, dim, delta) = (
        caps.next().unwrap()[0].parse::<u32>().unwrap(),
        caps.next().unwrap()[0].parse::<u64>().unwrap(),
        caps.next().unwrap()[0].parse::<u64>().unwrap(),
        caps.next().unwrap()[0].parse::<u32>().unwrap(),
        match &caps.next().unwrap()[0] {
            "1." => "1.0",
            "0.98999995." => "0.99",
            _ => unreachable!(),
        },
    );
    (algo, all_dim, seed, cnt, dim, delta)
}
fn main() {
    println!("algo,all_dim,seed,gen_cnt,dim,delta,sec,len,cnt");
    const BASE_PATH: &str = "results";

    let mut dirs: Vec<_> = read_dir(BASE_PATH, |t| t.is_dir()).unwrap().collect();
    dirs.sort();

    for dir in &dirs {
        let mut files: Vec<_> = read_dir(Path::new(BASE_PATH).join(dir), |t| t.is_file())
            .unwrap()
            .filter(|f| f.starts_with("ndim"))
            .collect();
        files.sort();

        for file in &files {
            let path = Path::new(BASE_PATH).join(dir).join(file);
            let s = cat(&path).unwrap();
            let (sec, len, cnt) = parse_sec_len_cnt(&s);
            let (algo, all_dim, seed, gen_cnt, dim, delta) = parse_path_to_tuple_random(&path);
            println!(
                "{},{},{},{},{},{},{},{},{}",
                algo, all_dim, seed, gen_cnt, dim, delta, sec, len, cnt
            );
        }
    }
}
