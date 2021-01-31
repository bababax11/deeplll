mod deeplll;
mod parse;

use deeplll::deep_lll;
use parse::matrix_parse;

use ndarray::prelude::*;
use std::path::Path;

use std::time::{Instant};

macro_rules! measure {
  ($str: expr, $x: expr) => {
    {
      let start = Instant::now();
      let result = $x;
      let end = start.elapsed();
      println!("{}: {}.{:03} sec", $str, end.as_secs(), end.subsec_nanos() / 1_000_000);
      result
    }
  };
}

fn main() {
  let path_str = "matrices/svp/svpchallengedim40seed0.txt";
  let path = Path::new(path_str);
  let b = matrix_parse(path);

  let (new_b, v, mu, hist, cnt) =
    measure!(path_str, deep_lll(b.slice(s![0..30, 0..30]).to_owned(), 1.into(), true, 100));
  println!("b: {:?}\nv_norms: {:?}\n{:?}\n(hist.len, cnt): {:?}\n{:?}", new_b, v, mu, (hist.len(), cnt), hist);
}
