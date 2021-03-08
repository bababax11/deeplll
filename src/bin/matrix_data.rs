use deeplll::deeplll::mat_to_str;
use deeplll::gen_mat::gen_mat;

const NDIMS: [usize; 7] = [10, 15, 20, 25, 30, 35, 40];

fn main() {
    const CNT: u64 = 1000;
    println!(";ndim;seed;matrix");

    for &ndim in &NDIMS {
        for seed in 0..5 {
            let mat = gen_mat(ndim, seed, CNT);
            print!(";{};{};{}", ndim, seed, mat_to_str(mat.view()));
        }
    }
}