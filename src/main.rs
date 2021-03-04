mod deeplll;
mod gen_mat;
mod parse;

use deeplll::{deep_lll_width, lll, mu::Mu, pot_lll, s2_lll, LLLFn};
use gen_mat::gen_mat;
use parse::matrix_parse;

use ndarray::prelude::*;
use rug::Rational;
use std::fs::File;
use std::io::prelude::*;
use std::path::Path;

use std::time::Instant;

macro_rules! measure {
    ($str: expr, $result_f: expr, $x: expr) => {{
        let start = Instant::now();
        let result = $x;
        let end = start.elapsed();
        writeln!(
            $result_f,
            "{}: {}.{:03} sec",
            $str,
            end.as_secs(),
            end.subsec_nanos() / 1_000_000
        )?;
        result
    }};
}

fn experiment_random(ndim: usize, seed: u64, cnt: usize) {
    let b = gen_mat(ndim, seed, cnt);
    let path_str_base = format!("ndim{}seed{}cnt{}", ndim, seed, cnt);
    experiment_mat(b, &[ndim], &path_str_base);
}

fn experiment_svp(mat_path_str: &str) {
    let path = Path::new(mat_path_str);
    let b = matrix_parse(path);
    let path_str_base = mat_path_str
        .split("/")
        .last()
        .unwrap()
        .split(".")
        .next()
        .unwrap();
    experiment_mat(b, &NDIMS, path_str_base);
}

fn experiment_mat(b: Array2<Rational>, ndims: &[usize], path_str_base: &str) {
    for &ndim in ndims {
        for rat in &[Rational::from(1), Rational::from((99, 100))] {
            macro_rules! experiment {
                ($dir_name: expr, $f: expr) => {
                    let result_path_str = format!(
                        "results/{}/{}dim{}delta{}.txt",
                        $dir_name,
                        path_str_base,
                        ndim,
                        rat.to_f32()
                    );
                    eprintln!("\n{}", &result_path_str);
                    experiment_unit(
                        b.slice(s![0..ndim, 0..ndim]),
                        &result_path_str,
                        rat.to_owned(),
                        $f,
                    )
                    .unwrap();
                };
            }

            experiment!("deeplll05", deep_lll_width(5));
            experiment!("deeplll10", deep_lll_width(10));

            // experiment!("lll", lll);

            // experiment!("potlll", pot_lll);

            // experiment!("s2lll", s2_lll);
        }
    }
}

fn experiment_unit<T: std::fmt::Debug>(
    b: ArrayView2<Rational>,
    result_path_str: &str,
    delta: Rational,
    f: impl LLLFn<T>,
) -> std::io::Result<()> {
    let mut result_f = File::create(result_path_str)?;
    let (new_b, v, mu, hist, cnt) =
        measure!(result_path_str, result_f, f(b.to_owned(), delta, true, 100));
    writeln!(
        result_f,
        "b: {:?}\nv_norms: {:?}\n{:?}\n(hist.len, cnt): {:?}\n{:?}",
        new_b,
        v,
        mu,
        (hist.len(), cnt),
        hist
    )?;
    Ok(())
}

const NDIMS: [usize; _] = [10, 15, 20, 25, 30, 35, 40];

fn main() {
    for i in 0..5 {
        let mat_path_str = format!("matrices/svp/svpchallengedim40seed{}.txt", i);
        experiment_svp(&mat_path_str);
    }
    const CNT: usize = 1000;

    for &ndim in &NDIMS {
      for seed in 0..5 {
        experiment_random(ndim, seed, CNT);
      }
    }
}
