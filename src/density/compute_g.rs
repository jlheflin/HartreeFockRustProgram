use ndarray::{Array2, Array4};


pub fn compute_g(dens_mat: &Array2<f64>, vee: &Array4<f64>) -> Array2<f64> {
    let nbasis = dens_mat.shape()[0];
    let mut g = Array2::<f64>::zeros((nbasis, nbasis));

    for i in 0..nbasis {
        for j in 0..nbasis {
            for k in 0..nbasis {
                for l in 0..nbasis {
                    let density = dens_mat[(k, l)];
                    let vee_ijkl = vee[(i, j, k, l)];
                    let vee_ilkj = vee[(i, l, k, j)];
                    g[(i, j)] += density * (vee_ijkl - 0.5 * vee_ilkj);
                }
            }
        }
    }
    return g;
}
