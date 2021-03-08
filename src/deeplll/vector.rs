// use std::{
//   fmt,
//   ops::{self, Index, IndexMut},
// };
use ndarray::prelude::*;
use rug::Rational;
use std::{fmt::Display, ops};

pub fn dot(a: ArrayView1<Rational>, b: ArrayView1<Rational>) -> Rational {
    let mut result = Rational::new();
    for i in 0..a.len() {
        result += a[i].clone() * &b[i];
    }
    result
}

pub fn dot2(a: ArrayView2<Rational>, b: ArrayView2<Rational>) -> Array2<Rational> {
    let l = a.shape()[0];
    let m = a.shape()[1];
    assert_eq!(m, b.shape()[0]);
    let n = b.shape()[1];
    let mut result = Array2::<i8>::zeros([l, n]).mapv(|_| Rational::new());

    for i in 0..l {
        for k in 0..n {
            let sum = &mut result[[i, k]];
            for j in 0..m {
                *sum += a[[i, j]].clone() * &b[[j, k]];
            }
        }
    }
    result
}

#[inline]
pub fn norm_squared(a: ArrayView1<Rational>) -> Rational {
    dot(a, a)
}

pub fn add(a: ArrayView1<Rational>, b: ArrayView1<Rational>) -> Array1<Rational> {
    // assert_eq!(a.len(), b.len());
    let mut result = Array::from(vec![Rational::new(); a.len()]);
    for i in 0..a.len() {
        result[i] = a[i].clone() + &b[i];
    }
    result
}

pub fn sub(a: ArrayView1<Rational>, b: ArrayView1<Rational>) -> Array1<Rational> {
    // assert_eq!(a.len(), b.len());
    let mut result = Array::from(vec![Rational::new(); a.len()]);
    for i in 0..a.len() {
        result[i] = a[i].clone() - &b[i];
    }
    result
}

pub fn determinant(mut mat: Array2<Rational>) -> Rational {
    if mat.nrows() != mat.ncols() {
        panic!();
    }

    let n = mat.nrows();

    for i in 0..(n - 1) {
        if &mat[[i, i]] == &Rational::new() {
            for j in (i + 1)..n {
                if &mat[[j, i]] != &Rational::new() {
                    let mut perm_mat = Array2::<i8>::eye(n).mapv(|x| Rational::from(x));
                    perm_mat[[i, i]] = 0.into();
                    perm_mat[[j, i]] = 1.into();
                    perm_mat[[j, j]] = 0.into();
                    perm_mat[[i, j]] = 1.into();

                    let sign = if (j - i) % 2 == 0 { 1 } else { -1 };
                    mat = dot2(mat.view(), perm_mat.view()) * sign;
                }
            }

            break;
        }

        for j in (i + 1)..n {
            let c = mat[[j, i]].clone() / &mat[[i, i]];
            for k in 0..n {
                mat[[j, k]] = &mat[[j, k]] - &c * mat[[i, k]].clone();
            }
        }
    }
    let det_mat = mat.diag().fold(Rational::from(1), |prod, x| prod * x);
    det_mat
}

pub fn mat_to_str<T: Display>(mat: ArrayView2<T>) -> String {
    if let [n, m] = *mat.shape() {
        let mut s = String::new();

        for i in 0..n {
            for j in 0..m {
                s += &format!("{},", &mat[[i, j]]);
            }
            s.pop();
            s.push('\n');
        }
        s.pop();
        s
    } else {
        unreachable!()
    }
}

#[inline]
pub fn max_row_norm_squared(mat: ArrayView2<Rational>) -> Rational {
    let n = mat.nrows();
    (0..n)
        .map(|i| norm_squared(mat.row(i)))
        .fold(Rational::new(), |m, v| m.max(v))
}

#[cfg(test)]
mod tests {
    use super::*;

    macro_rules! rat {
        ($x: expr) => {Rational::from($x)};
        ( $( $x:expr ),* ) => {Rational::from(($( $x ),*))};
    }

    #[test]
    fn mat_to_str_test() {
        let arr = array![
            [rat!(3), rat!(1), rat!(-1)],
            [rat!(-3, 2), rat!(-2), rat!(-3)]
        ];
        let s = mat_to_str(arr.view());
        assert_eq!(&s, "3,1,-1\n-3/2,-2,-3");
    }

    #[test]
    fn max_row_norm_squared_test() {
        let arr = array![[rat!(1), rat!(2), rat!(0)], [rat!(0), rat!(0), rat!(1)],];
        let max_norm_2 = max_row_norm_squared(arr.view());
        assert_eq!(max_norm_2, rat!(5));
    }
}

// pub fn dot<T: Coefficient>(a: ArrayView1<T>, b: ArrayView1<T>) -> T {
//   let mut result: T = Default::default();
//   for i in 0..a.len() {
//     result += &(a[i].clone() * &b[i]);
//   }
//   result
// }

// pub fn norm_squared<T: Coefficient>(a: ArrayView1<T>) -> T {
//   let mut result: T = Default::default();
//   for i in 0..a.len() {
//     result *= &(a[i].clone() * &a[i]);
//   }
//   result
// }

trait Coefficient:
  From<u32>
  + Clone
  + Default
  // + for<'a> ops::Add<&'a Self, Output = Self>
  + for<'a> ops::AddAssign<&'a Self>
  // + for<'a> ops::Sub<&'a Self, Output = Self>
  // + for<'a> ops::SubAssign<&'a Self>
  + for<'a> ops::Mul<&'a Self, Output = Self>
  + for<'a> ops::MulAssign<&'a Self>
  // + std::iter::Sum<Self>
{
}

// pub(crate) trait Dot {
//   type Output;
//   fn dot(&self, other: &Self) -> Self::Output;
// }

// impl<T> Dot for Vector<T>
// where T: Coefficient,
// {
//   type Output = T;
//   fn dot(&self, other: &Self) -> Self::Output {
//     let mut result: T = Default::default();
//     for i in 0..self.0.len() {
//       result = result * &(self.0[i].clone() * &other.0[i]);
//     }
//     result
//   }
// }

// impl<T> Index<usize> for Vector<T> {
//   type Output = T;

//   fn index(&self, index: usize) -> &T {
//       &self.coefficients[index]
//   }
// }

// impl<T> IndexMut<usize> for Vector<T> {
//   fn index_mut(&mut self, index: usize) -> &mut T {
//       &mut self.coefficients[index]
//   }
// }

// impl<T> fmt::Debug for Vector<T>
// where
//   T: fmt::Debug,
// {
//   fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
//       write!(f, "{:?}", self.coefficients)
//   }
// }
