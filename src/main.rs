mod deeplll;
mod parse;

use deeplll::deep_lll;
use parse::matrix_parse;

use ndarray::prelude::*;
use std::path::Path;

fn main() {
  let path = Path::new("matrices/svp/svpchallengedim40seed0.txt");
  let b = matrix_parse(&path);
  let (new_b, v, mu, hist, cnt) =
    deep_lll(b.slice(s![0..10, 0..10]).to_owned(), 1.into(), true, 100);
  dbg!((new_b, v, mu, hist, cnt));
}
