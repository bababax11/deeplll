mod deeplll;
mod parse;

use deeplll::{deep_lll, s2_lll};
use parse::matrix_parse;

use ndarray::prelude::*;
use rug::Rational;
use std::path::Path;
use std::fs::File;
use std::io::prelude::*;

use std::time::{Instant};

macro_rules! measure {
  ($str: expr, $result_f: expr, $x: expr) => {
    {
      let start = Instant::now();
      let result = $x;
      let end = start.elapsed();
      write!($result_f, "{}: {}.{:03} sec", $str, end.as_secs(), end.subsec_nanos() / 1_000_000)?;
      result
    }
  };
}

fn experiment(b: Array2<Rational>, result_path_str: &str) -> std::io::Result<()> {
  let mut result_f = File::create(result_path_str)?;
  let (new_b, v, mu, hist, cnt) =
    measure!(result_path_str, result_f, deep_lll(b.slice(s![0..10, 0..10]).to_owned(), 1.into(), true, 100));
  writeln!(result_f, "b: {:?}\nv_norms: {:?}\n{:?}\n(hist.len, cnt): {:?}\n{:?}", new_b, v, mu, (hist.len(), cnt), hist)?;
  Ok(())
}

fn main() {
  let path_str = "matrices/svp/svpchallengedim40seed0.txt";
  let path = Path::new(path_str);
  let b = matrix_parse(path);

  experiment(b, "results/svp/svpchallengedim40seed0.txt").unwrap();
}
