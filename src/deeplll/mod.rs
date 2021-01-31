mod vector;
use vector::{dot, norm_squared};
// mod matrix;
// use matrix::Matrix;
mod mu;
use mu::Mu;
use rug::{Rational, ops::Pow};
use ndarray::prelude::*;

fn make_v(tmp_v: Array2<Rational>) -> Array1<Rational> {
  let n = tmp_v.nrows();
  let mut v = Array::from(
    vec![Rational::default(); n]
  );
  for i in 0..n {
    v[i] = dot(tmp_v.row(i), tmp_v.row(i));
  }
  v
}

fn orthogonize(b: &Array2<Rational>) -> (Array1<Rational>, Mu) {
  let n = b.nrows();
  let mut tmp_v = b.clone();
  let mut mu = Mu::new(n);
  for i in 0..n {
    for j in 0..i {
      let tmp = dot(b.row(i), tmp_v.row(j)) / norm_squared(tmp_v.row(j));
      let s = vector::sub(tmp_v.row(i),
       tmp_v.row(j).map(|x| Rational::from(&tmp * x)).view());
      tmp_v.row_mut(i).assign(&s);
      mu[(i, j)] = tmp;
    }
  }
  let v = make_v(tmp_v);
  (v, mu)
}

fn size_reduce(b: &Array2<Rational>, mu: &Mu) -> (Array2<Rational>, Mu) {
  let n = b.nrows();
  let mut new_b = b.clone();
  for k in 0..n {
    for j in 0..k {
      let s  = vector::sub(new_b.row(k),
      b.row(j).map(|x| x * Rational::from(mu[(k, j)].round_ref())).view());
      new_b.row_mut(k).assign(&s);
    }
  }
  let (_, new_mu) = orthogonize(&new_b);
  (new_b, new_mu)
}

fn swap(mut b: ArrayViewMut2<Rational>, i: usize, k: usize) {
  let row_i = b.slice(s![i, ..]).to_owned();
  let tmp = b.slice(s![k, ..]).to_owned();
  b.slice_mut(s![i, ..]).assign(&tmp);
  for l in ((i+1)..k).rev() {
    let tmp = b.slice(s![l, ..]).to_owned();
    b.slice_mut(s![l + 1, ..]).assign(&tmp);
  }
  b.slice_mut(s![i + 1, ..]).assign(&row_i);
}

pub fn deep_lll(mut b: Array2<Rational>, delta: Rational, verbose: bool, verbose_count: usize) -> (Array2<Rational>, Array1<Rational>, Mu, Vec<(usize, usize)>, usize) {
  let n = b.nrows();
  let mut hist = Vec::with_capacity(verbose_count * 5);
  let mut k = 1;
  let (mut v, mut mu) = orthogonize(&b);

  let mut cnt = 0;
  while k < n {
    let (_b, _mu) = size_reduce(&b, &mu);
    b = _b;
    mu = _mu;

    let mut after_break = false;
    let mut c = norm_squared(b.row(k));

    for i in 0..k {

      cnt += 1;
      if verbose && cnt % verbose_count == 0 {
        print!("=");
      }

      if c >= Rational::from(&delta * &v[i]) {
        c -= Rational::from((&mu[(k, i)]).pow(2)) * &v[i];
      } else {
        hist.push((i, k));
        swap(b.view_mut(), i, k);
        let (_v, _mu) = orthogonize(&b);
        v = _v;
        mu = _mu;
        k = usize::max(i, 1);
        after_break = true;
        break;
      }
    }
    if !after_break {
      k += 1;
    }
  }
  (b, v, mu, hist, cnt)
}
#[cfg(test)]
mod tests {
  use super::*;

  macro_rules! rat {
    ($x: expr) => {Rational::from($x)};
    ( $( $x:expr ),* ) => {Rational::from(($( $x ),*))};
  }
  #[test]
  fn orthogonize_test() {
    let b = array![
      [rat!(1), rat!(1), rat!(0)],
      [rat!(0), rat!(1), rat!(1)],
      [rat!(1), rat!(0), rat!(1)]
    ];
    let (v, m) = orthogonize(&b);
    assert_eq!(v, array![
      rat!(2), rat!(3, 2), rat!(4, 3)
    ]);
    assert_eq!(m[(1, 0)], rat!(1, 2));
    assert_eq!(m[(2, 0)], rat!(1, 2));
    assert_eq!(m[(2, 1)], rat!(1, 3));
  }
  #[test]
  fn size_reduce_test() {
    let b = array![
      [rat!(3), rat!(1), rat!(-1)],
      [rat!(-3), rat!(-2), rat!(-3)],
      [rat!(1), rat!(-1), rat!(2)]
    ];
    let (_, mu) = orthogonize(&b);
    assert_eq!(mu[(1, 0)], rat!(-8, 11));
    assert_eq!(mu[(2, 0)], rat!(0));
    assert_eq!(mu[(2, 1)], rat!(-77, 178));
    let (b2, mu2) = size_reduce(&b, &mu);
    assert_eq!(b2, array![
      [rat!(3), rat!(1), rat!(-1)],
      [rat!(0), rat!(-1), rat!(-4)],
      [rat!(1), rat!(-1), rat!(2)]
    ]);
    assert_eq!(mu2[(1, 0)], rat!(3, 11));
    assert_eq!(mu2[(2, 0)], rat!(0));
    assert_eq!(mu2[(2, 1)], rat!(-77, 178));
  }
  #[test]
  fn deep_lll_test() {
    let b = array![
      [rat!(0), rat!(3), rat!(-2)],
      [rat!(-3), rat!(1), rat!(-2)],
      [rat!(2), rat!(-2), rat!(-2)],
    ];
    let (b2, v2, mu2, hist, _) = deep_lll(b, rat!(1), false, 10);
    assert_eq!(b2, array![
      [rat!(2), rat!(-2), rat!(-2)],
      [rat!(0), rat!(3), rat!(-2)],
      [rat!(-3), rat!(-2), rat!(0)]
    ]);
    assert_eq!(v2, array![
      rat!(12), rat!(38, 3), rat!(19, 2)
    ]);
    assert_eq!(mu2[(1, 0)], rat!(-1, 6));
    assert_eq!(mu2[(2, 0)], rat!(-1, 6));
    assert_eq!(mu2[(2, 1)], rat!(-1, 2));
    assert_eq!(&hist, &[(0, 2)]);
  }
}