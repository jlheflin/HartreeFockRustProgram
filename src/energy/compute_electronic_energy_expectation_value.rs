use ndarray::Array2;


pub fn compute_electronic_energy_expectation_value(dens_mat: &Array2<f64>, t: &Array2<f64>, vne: &Array2<f64>, g: &Array2<f64>) -> f64 {
    let mut electronic_energy = 0.0;
    
    let hcore = t + vne;
    let nbasis = dens_mat.shape()[0];

    for i in 0..nbasis {
        for j in 0..nbasis {
            electronic_energy += dens_mat[(i, j)] * (hcore[(i, j)] + 0.5 * g[(i, j)]);
        }
    }
    return electronic_energy;
}
