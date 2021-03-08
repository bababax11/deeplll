use deeplll::gen_mat::gen_mat;
use deeplll::{
    deeplll::vector::{determinant, max_row_norm_squared},
    parse::matrix_parse,
};
use ndarray::s;

const NDIMS: [usize; 7] = [10, 15, 20, 25, 30, 35, 40];

fn main() {
    const CNT: u64 = 1000;
    println!(
        "type,ndim,seed,max_row_norm_squared,det,max_row_norm_squared_f64,max_row_norm_f64,det_f64"
    );

    for seed in 0..5 {
        let mat_path_str = format!("matrices/svp/svpchallengedim40seed{}.txt", seed);
        let mat = matrix_parse(mat_path_str.as_ref());

        for &ndim in &NDIMS {
            let part_mat = mat.slice(s![0..ndim, 0..ndim]);
            let max_norm_2 = max_row_norm_squared(part_mat);
            let det = determinant(part_mat.to_owned());
            println!(
                "{},{},{},{},{},{:e},{:e},{:e}",
                "svp",
                ndim,
                seed,
                max_norm_2,
                det,
                max_norm_2.to_f64(),
                max_norm_2.to_f64().sqrt(),
                det.to_f64()
            );
        }
    }

    for &ndim in &NDIMS {
        for seed in 0..5 {
            let mat = gen_mat(ndim, seed, CNT);
            let max_norm_2 = max_row_norm_squared(mat.view());
            let det = determinant(mat);
            println!(
                "{},{},{},{},{},{:e},{:e},{:e}",
                "random",
                ndim,
                seed,
                max_norm_2,
                det,
                max_norm_2.to_f64(),
                max_norm_2.to_f64().sqrt(),
                det.to_f64()
            );
        }
    }
}
