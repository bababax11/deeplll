use std::path::Path;

use deeplll::parse::{cat, parse_sec_len_cnt, read_dir};

fn main() {
    println!("file_path,sec,len,cnt");
    const BASE_PATH: &str = "results";

    let mut dirs: Vec<_> = read_dir(BASE_PATH, |t| t.is_dir()).unwrap().collect();
    dirs.sort();

    for dir in &dirs {
        let mut files: Vec<_> = read_dir(Path::new(BASE_PATH).join(dir), |t| t.is_file())
            .unwrap()
            .collect();
        files.sort();

        for file in &files {
            let path = Path::new(BASE_PATH).join(dir).join(file);
            let s = cat(&path).unwrap();
            let (sec, len, cnt) = parse_sec_len_cnt(&s);
            println!("{},{},{},{}", path.to_str().unwrap(), sec, len, cnt);
        }
    }
}
