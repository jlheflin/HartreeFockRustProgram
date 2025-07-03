use ndarray::Array2;

pub fn compute_density_matrix(mos: Array2<f64>) -> Array2<f64> {
    let nbasis = mos.shape()[0];
    let mut dens_mat = Array2::<f64>::zeros((nbasis, nbasis));

    let occupation = 2;
    let number_oo = 1;

    for i in 0..nbasis {
        for j in 0..nbasis {
            for oo in 0..number_oo {
                let c = mos[(i, oo)];
                let c_dagger = mos[(j, oo)];
                dens_mat[(i, j)] += occupation as f64 * c * c_dagger;
            }
        }
    }
    return dens_mat;
}
