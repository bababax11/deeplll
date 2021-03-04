use ndarray::prelude::*;
use regex::Regex;
use rug::Rational;
use std::fs::File;
use std::io;
use std::io::prelude::*;
use std::path::Path;

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

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn matrix_parse_test() {
        let path = Path::new("matrices/svp/svpchallengedim40seed0.txt");
        let mat = matrix_parse(path);
        dbg!(mat);
    }
}
