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
      writeln!($result_f, "{}: {}.{:03} sec", $str, end.as_secs(), end.subsec_nanos() / 1_000_000)?;
      result
    }
  };
}

fn experiment_mat(mat_path_str: &str, dir_name: &str) {
  assert_eq!(dir_name, "deeplll"); // TODO!!!!!!!!!!!!!!!!!!
  let path = Path::new(mat_path_str);
  let b = matrix_parse(path);
  for ndim in &[15, 20, 25, 30] {
    let path_str_base = mat_path_str.split("/").last().unwrap().split(".").next().unwrap();
    for rat in &[Rational::from(1), Rational::from((99, 100))] {
      let result_path_str = format!("results/{}/{}dim{}delta{}.txt", dir_name, path_str_base, *ndim, rat.to_f32());
      eprintln!("\n{}", &result_path_str);
      experiment_unit_deeplll(b.slice(s![0..(*ndim), 0..(*ndim)]), &result_path_str, rat.to_owned()).unwrap();
    }
  }
}

fn experiment_unit_deeplll(b: ArrayView2<Rational>, result_path_str: &str, delta: Rational) -> std::io::Result<()> {
  let mut result_f = File::create(result_path_str)?;
  let (new_b, v, mu, hist, cnt) =
    measure!(result_path_str, result_f, deep_lll(b.to_owned(), delta, true, 100));
  writeln!(result_f, "b: {:?}\nv_norms: {:?}\n{:?}\n(hist.len, cnt): {:?}\n{:?}", new_b, v, mu, (hist.len(), cnt), hist)?;
  Ok(())
}

fn main() {
  const DIR_NAME: &str = "deeplll";
  for i in 0..5 {
    let mat_path_str = format!("matrices/svp/svpchallengedim40seed{}.txt", i);
    experiment_mat(&mat_path_str, DIR_NAME);
  }
}
