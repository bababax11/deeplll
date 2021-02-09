mod vector;
use vector::{dot, norm_squared};
// mod matrix;
// use matrix::Matrix;
pub mod mu;
use mu::Mu;
use ndarray::prelude::*;
use rug::{ops::Pow, Rational};

fn make_v(tmp_v: Array2<Rational>) -> Array1<Rational> {
  let n = tmp_v.nrows();
  let mut v = Array::from(vec![Rational::default(); n]);
  for i in 0..n {
    v[i] = dot(tmp_v.row(i), tmp_v.row(i));
  }
  v
}

fn gso(b: &Array2<Rational>) -> (Array1<Rational>, Mu) {
  let n = b.nrows();
  let mut tmp_v = b.clone();
  let mut mu = Mu::new(n);
  for i in 0..n {
    for j in 0..i {
      let tmp = dot(b.row(i), tmp_v.row(j)) / norm_squared(tmp_v.row(j));
      let s = vector::sub(
        tmp_v.row(i),
        tmp_v.row(j).map(|x| &tmp * x.clone()).view(),
      );
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
      let s = vector::sub(
        new_b.row(k),
        b.row(j)
          .map(|x| x * (Rational::from((-1, 2)) + &mu[(k, j)]).ceil())
          .view(),
      );
      new_b.row_mut(k).assign(&s);
    }
  }
  let (_, new_mu) = gso(&new_b);
  (new_b, new_mu)
}

fn swap(mut b: ArrayViewMut2<Rational>, i: usize, k: usize) {
  let row_i = b.slice(s![i, ..]).to_owned();
  let tmp = b.slice(s![k, ..]).to_owned();
  b.slice_mut(s![i, ..]).assign(&tmp);
  for l in ((i + 1)..k).rev() {
    let tmp = b.slice(s![l, ..]).to_owned();
    b.slice_mut(s![l + 1, ..]).assign(&tmp);
  }
  b.slice_mut(s![i + 1, ..]).assign(&row_i);
}

pub fn lll(
  mut b: Array2<Rational>,
  delta: Rational,
  verbose: bool,
  verbose_count: usize,
) -> (
  Array2<Rational>,
  Array1<Rational>,
  Mu,
  Vec<(usize, usize)>,
  usize,
) {
  let n = b.nrows();
  let mut hist = Vec::with_capacity(verbose_count * 100);
  let mut k = 1;
  let (mut v, mut mu) = gso(&b);
  
  let mut cnt = 0;
  while k < n {
    cnt += 1;
    if verbose && cnt % verbose_count == 0 {
      eprint!("=");
    }
    let (_b, _mu) = size_reduce(&b, &mu);
    b = _b;
    mu = _mu;
    
    if &v[k] >= &((&delta - mu[(k, k-1)].clone().pow(2)) * &v[k-1]) {
      k += 1;
    } else {
      hist.push((k - 1, k));
      swap(b.view_mut(), k - 1, k);
      let (_v, _mu) = gso(&b);
      v = _v;
      mu = _mu;
      k = usize::max(k - 1, 1);
    }
  }
  (b, v, mu, hist, cnt)
}

pub fn deep_lll(
  mut b: Array2<Rational>,
  delta: Rational,
  verbose: bool,
  verbose_count: usize,
) -> (
  Array2<Rational>,
  Array1<Rational>,
  Mu,
  Vec<(usize, usize)>,
  usize,
) {
  let n = b.nrows();
  let mut hist = Vec::with_capacity(verbose_count * 100);
  let mut k = 1;
  let (mut v, mut mu) = gso(&b);

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
        eprint!("=");
      }

      if c >= delta.clone() * &v[i] {
        c -= mu[(k, i)].clone().pow(2) * &v[i];
      } else {
        hist.push((i, k));
        swap(b.view_mut(), i, k);
        let (_v, _mu) = gso(&b);
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

fn max_sjk(v: &Array1<Rational>, mu: &Mu, k: usize) -> (usize, Rational) {
  let mut sjk = vec![Rational::new(); k];
  let mut before = Rational::new();
  let mut pjbk = v[k].clone();
  for (e, j) in (0..k).rev().enumerate() {
    let tmp = mu[(k, j)].clone().pow(2) * &v[j];
    pjbk += &tmp;
    before += tmp * (v[j].clone() / &pjbk - 1);
    sjk[e] = before.clone();
  }
  let (e_, max_val) = argmax(&sjk);
  (k - e_ - 1, max_val)
}
fn min_pot(v: &Array1<Rational>, mu: &Mu, k: usize) -> (usize, Rational) {
  let mut potjk = vec![Rational::new(); k];
  let mut before = Rational::from(1);
  let mut pjbk = v[k].clone();
  for (e, j) in (0..k).rev().enumerate() {
    let tmp = mu[(k, j)].clone().pow(2) * &v[j];
    pjbk += &tmp;
    before *= pjbk.clone() / &v[j];
    potjk[e] = before.clone();
  }
  let (e_, min_val) = argmin(&potjk);
  (k - e_ - 1, min_val)
}

fn argmax<T: PartialOrd + Clone>(xs: &[T]) -> (usize, T) {
  if xs.is_empty() {
    panic!();
  } else if xs.len() == 1 {
    (0, xs[0].clone())
  } else {
    let mut max_val = &xs[0];
    let mut max_ix = 0;
    for (i, x) in xs.iter().enumerate().skip(1) {
      if x > max_val {
        max_val = x;
        max_ix = i;
      }
    }
    (max_ix, max_val.clone())
  }
}
fn argmin<T: PartialOrd + Clone>(xs: &[T]) -> (usize, T) {
  if xs.is_empty() {
    panic!();
  } else if xs.len() == 1 {
    (0, xs[0].clone())
  } else {
    let mut min_val = &xs[0];
    let mut min_ix = 0;
    for (i, x) in xs.iter().enumerate().skip(1) {
      if x < min_val {
        min_val = x;
        min_ix = i;
      }
    }
    (min_ix, min_val.clone())
  }
}
#[inline]
pub fn ss(v: &Array1<Rational>) -> Rational {
  v.fold(Rational::new(), |sum, e| sum + e)
}
pub fn pot(v: &Array1<Rational>) -> Rational {
  let mut p = Rational::from(1);
  let n = v.len();
  for (i, r) in v.iter().enumerate() {
    p *= Rational::from(r.pow((n - i) as u32));
  }
  p
}

pub fn s2_lll(
  mut b: Array2<Rational>,
  delta: Rational,
  verbose: bool,
  verbose_count: usize,
) -> (
  Array2<Rational>,
  Array1<Rational>,
  Mu,
  Vec<(usize, usize, Rational)>,
  usize,
) {
  let n = b.nrows();
  let mut hist = Vec::with_capacity(verbose_count * 100);
  let mut k = 1;
  let (mut v, mut mu) = gso(&b);

  let mut cnt = 0;

  while k < n {
    cnt += 1;
    if verbose && cnt % verbose_count == 0 {
      eprint!("=");
    }
    let (_b, _mu) = size_reduce(&b, &mu);
    b = _b;
    mu = _mu;

    let (i, s) = max_sjk(&v, &mu, k);
    if s <= (Rational::from(1) - &delta) * ss(&v) {
      k += 1;
    } else {
      hist.push((i, k, ss(&v)));
      swap(b.view_mut(), i, k);
      let (_v, _mu) = gso(&b);
      v = _v;
      mu = _mu;
      k = usize::max(i, 1);
    }
  }
  (b, v, mu, hist, cnt)
}

pub fn pot_lll(
  mut b: Array2<Rational>,
  delta: Rational,
  verbose: bool,
  verbose_count: usize,
) -> (
  Array2<Rational>,
  Array1<Rational>,
  Mu,
  Vec<(usize, usize)>,
  usize,
) {
  let n = b.nrows();
  let mut hist = Vec::with_capacity(verbose_count * 100);
  let mut k = 1;
  let (mut v, mut mu) = gso(&b);

  let mut cnt = 0;

  while k < n {
    cnt += 1;
    if verbose && cnt % verbose_count == 0 {
      eprint!("=");
    }
    let (_b, _mu) = size_reduce(&b, &mu);
    b = _b;
    mu = _mu;

    let (i, p) = min_pot(&v, &mu, k);
    if &delta <= &p {
      k += 1;
    } else {
      hist.push((i, k));
      swap(b.view_mut(), i, k);
      let (_v, _mu) = gso(&b);
      v = _v;
      mu = _mu;
      k = usize::max(i, 1);
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
  fn gso_test() {
    let b = array![
      [rat!(1), rat!(1), rat!(0)],
      [rat!(0), rat!(1), rat!(1)],
      [rat!(1), rat!(0), rat!(1)]
    ];
    let (v, m) = gso(&b);
    assert_eq!(pot(&v), rat!(24));
    assert_eq!(v, array![rat!(2), rat!(3, 2), rat!(4, 3)]);
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
    let (_, mu) = gso(&b);
    assert_eq!(mu[(1, 0)], rat!(-8, 11));
    assert_eq!(mu[(2, 0)], rat!(0));
    assert_eq!(mu[(2, 1)], rat!(-77, 178));
    let (b2, mu2) = size_reduce(&b, &mu);
    assert_eq!(
      b2,
      array![
        [rat!(3), rat!(1), rat!(-1)],
        [rat!(0), rat!(-1), rat!(-4)],
        [rat!(1), rat!(-1), rat!(2)]
      ]
    );
    assert_eq!(mu2[(1, 0)], rat!(3, 11));
    assert_eq!(mu2[(2, 0)], rat!(0));
    assert_eq!(mu2[(2, 1)], rat!(-77, 178));
  }
  #[test]
  fn lll_test() {
    let b = array![
      [rat!(1), rat!(1), rat!(1)],
      [rat!(-1), rat!(0), rat!(2)],
      [rat!(3), rat!(5), rat!(6)]
    ];
    let (b2, _, _, _, _) = lll(b.clone(), rat!(1), false, 0);
    assert_eq!(b2, array![
      [rat!(0), rat!(1), rat!(0)],
      [rat!(1), rat!(0), rat!(1)],
      [rat!(-1), rat!(0), rat!(2)]
    ]);
    let (b2, _, _, _, _) = lll(b, rat!(75, 100), false, 0);
    assert_eq!(b2, array![
      [rat!(0), rat!(1), rat!(0)],
      [rat!(1), rat!(0), rat!(1)],
      [rat!(-1), rat!(0), rat!(2)]
    ]);
  }
  #[test]
  fn deep_lll_test() {
    let b = array![
      [rat!(0), rat!(3), rat!(-2)],
      [rat!(-3), rat!(1), rat!(-2)],
      [rat!(2), rat!(-2), rat!(-2)],
    ];
    let (b2, v2, mu2, hist, _) = deep_lll(b, rat!(1), false, 0);
    assert_eq!(
      b2,
      array![
        [rat!(2), rat!(-2), rat!(-2)],
        [rat!(0), rat!(3), rat!(-2)],
        [rat!(-3), rat!(1), rat!(-2)]
      ]
    );
    assert_eq!(v2, array![rat!(12), rat!(38, 3), rat!(19, 2)]);
    assert_eq!(mu2[(1, 0)], rat!(-1, 6));
    assert_eq!(mu2[(2, 0)], rat!(-1, 3));
    assert_eq!(mu2[(2, 1)], rat!(1, 2));
    assert_eq!(&hist, &[(0, 2)]);
  }
  #[test]
  fn s2_lll_test() {
    let b = array![
      [rat!(3), rat!(1), rat!(-1)],
      [rat!(0), rat!(-1), rat!(-4)],
      [rat!(1), rat!(-1), rat!(2)]
    ];
    let (b2, v2, mu2, hist, _) = s2_lll(b, rat!(1), false, 0);
    assert_eq!(b2, array![
      [rat!(1), rat!(-1), rat!(2)],
      [rat!(1), rat!(-2), rat!(-2)],
      [rat!(3), rat!(1), rat!(-1)]
    ]);
    assert_eq!(v2, array![
      rat!(6), rat!(53, 6), rat!(529, 53) 
    ]);
    assert_eq!(mu2[(1, 0)], rat!(-1, 6));
    assert_eq!(mu2[(2, 0)], rat!(0));
    assert_eq!(mu2[(2, 1)], rat!(18, 53));
    assert_eq!(&hist, &[(1, 2, rat!(59041, 1958)), (0, 2, rat!(1651, 66)), (0, 2, rat!(2239, 90))]);
  }
  #[test]
  fn pot_lll_test() {
    let b = array![
      [rat!(3), rat!(1), rat!(-1)],
      [rat!(-3), rat!(-2), rat!(-3)],
      [rat!(1), rat!(-1), rat!(2)]
    ];
    let (b2, v2, mu2, _, _) = pot_lll(b, rat!(1), false, 0);
    assert_eq!(b2, array![
      [rat!(1), rat!(-1), rat!(2)],
      [rat!(1), rat!(-2), rat!(-2)],
      [rat!(3), rat!(1), rat!(-1)]
    ]);
    assert_eq!(v2, array![
      rat!(6), rat!(53, 6), rat!(529, 53)
    ]);
    assert_eq!(mu2[(1, 0)], rat!(-1, 6));
    assert_eq!(mu2[(1, 1)], rat!(0));
    assert_eq!(mu2[(1, 2)], rat!(18, 53));
  }
}
