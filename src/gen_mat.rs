use crate::deeplll::vector::add;
use ndarray::prelude::*;
use rand::Rng;
use rug::Rational;

pub fn gen_mat(ndim: usize, seed: u64, cnt: u64) -> Array2<Rational> {
    use rand::SeedableRng;
    let mut rng = rand_xoshiro::Xoshiro256StarStar::seed_from_u64(seed);

    let mut arr = Array2::<i8>::eye(ndim).mapv(|x| Rational::from(x));
    for _ in 0..cnt {
        let mut i;
        let mut j;
        loop {
            i = rng.gen_range(0..ndim);
            j = rng.gen_range(0..ndim);
            if i != j {
                break;
            }
        }
        let tmp = add(
            arr.slice(s![i, ..]),
            (arr.slice(s![j, ..]).to_owned() * (-1i32).pow(rng.gen_range(0..2))).view(),
        );
        arr.slice_mut(s![i, ..]).assign(&tmp);
    }
    arr
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::deeplll::vector::determinant;

    #[test]
    fn gen_mat_test() {
        for i in 0..3 {
            let mat = gen_mat(10, i, 1000);
            dbg!(&mat);
            assert_eq!(determinant(mat), 1);
        }
    }
}
