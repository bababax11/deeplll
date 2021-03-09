use ndarray::prelude::*;
use once_cell::sync::Lazy;
use regex::Regex;
use rug::Rational;
use std::fs::File;
use std::io::prelude::*;
use std::path::Path;
use std::{fs, io};

// A simple implementation of `% cat path`
// `% cat path`のシンプルな実装
pub fn cat(path: &Path) -> io::Result<String> {
    let mut f = File::open(path)?;
    let mut s = String::new();
    match f.read_to_string(&mut s) {
        Ok(_) => Ok(s),
        Err(e) => Err(e),
    }
}

pub fn read_dir<P: AsRef<Path>>(
    path: P,
    file_or_dir: impl Fn(&fs::FileType) -> bool,
) -> io::Result<impl Iterator<Item = String>> {
    Ok(fs::read_dir(path)?.filter_map(move |entry| {
        let entry = entry.ok()?;
        if file_or_dir(&entry.file_type().ok()?) {
            Some(entry.file_name().to_string_lossy().into_owned())
        } else {
            None
        }
    }))
}

pub fn matrix_parse(path: &Path) -> Array2<Rational> {
    let mut nrows = 0;
    let s = cat(&path).unwrap();
    let re_line = Regex::new(r"\[[\s\d]+\]").unwrap();
    let re_num = Regex::new(r"\d+").unwrap();
    // let caps = re.captures_iter(&s);
    let mut rat_vec = Vec::with_capacity(1600);
    for caps in re_line.captures_iter(&s) {
        nrows += 1;

        for num in re_num.captures_iter(&caps[0]) {
            let rat = num[0].parse::<Rational>().unwrap();
            rat_vec.push(rat);
        }
    }
    assert_eq!(rat_vec.len() % nrows, 0);
    Array::from_shape_vec((nrows, rat_vec.len() / nrows), rat_vec).unwrap()
}

pub fn parse_sec_len_cnt(s: &str) -> (f32, u32, u32) {
    static RE_SEC: Lazy<Regex> = Lazy::new(|| Regex::new(r"(.[0-9])+ sec").unwrap());
    static RE_LEN_CNT: Lazy<Regex> =
        Lazy::new(|| Regex::new(r"\(hist.len, cnt\): \([0-9]+, [0-9]+\)").unwrap());
    let sec_str = &RE_SEC.captures(s).unwrap()[0];
    let mut sec_splited = sec_str.split(' ');
    let sec = if let Ok(sec) = sec_splited.next().unwrap().parse::<f32>() {
        sec
    } else {
        sec_splited.next().unwrap().parse::<f32>().unwrap()
    };

    let len_cnt = &RE_LEN_CNT.captures(s).unwrap()[0];

    static RE_TUP: Lazy<Regex> = Lazy::new(|| Regex::new(r"[0-9]+").unwrap());
    let len_cnt = &mut RE_TUP.captures_iter(len_cnt);
    let (len, cnt) = (&len_cnt.next().unwrap()[0], &len_cnt.next().unwrap()[0]);
    (
        sec,
        len.parse::<u32>().unwrap(),
        cnt.parse::<u32>().unwrap(),
    )
}

#[cfg(test)]
mod tests {
    use average::assert_almost_eq;

    use super::*;
    #[test]
    fn matrix_parse_test() {
        let path = Path::new("matrices/svp/svpchallengedim40seed0.txt");
        let mat = matrix_parse(path);
        dbg!(mat);
    }
    #[test]
    fn parse_sec_len_cnt_test() {
        let s = cat("results/deeplll/ndim10seed0cnt1000dim10delta0.98999995.txt".as_ref()).unwrap();
        let (sec, len, cnt) = parse_sec_len_cnt(&s);

        assert_almost_eq!(sec, 3.951, 1e-5);
        assert_eq!(len, 486);
        assert_eq!(cnt, 6316);
    }
}
