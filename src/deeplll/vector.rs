// use std::{
//   fmt,
//   ops::{self, Index, IndexMut},
// };
use std::ops;
use rug::Rational;
use ndarray::prelude::*;

pub fn dot(a: ArrayView1<Rational>, b: ArrayView1<Rational>) -> Rational {
  let mut result = Rational::default();
  for i in 0..a.len() {
    result += &(a[i].clone() * &b[i]);
  }
  result
}

pub fn norm_squared(a: ArrayView1<Rational>) -> Rational {
  let mut result = Rational::default();
  for i in 0..a.len() {
    result += &(a[i].clone() * &a[i]);
  }
  result
}

pub fn sub(a: ArrayView1<Rational>, b: ArrayView1<Rational>) -> Array1<Rational> {
  assert_eq!(a.len(), b.len());
  let mut result = Array::from(
    vec![Rational::default(); a.len()]
  );
  for i in 0..a.len() {
    result[i] = a[i].clone() - &b[i];
  }
  result
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