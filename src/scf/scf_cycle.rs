use ndarray::{Array2, Array4};
use ndarray_linalg::{Eigh, UPLO};

use crate::{classes::molecule::Molecule, density::{compute_density_matrix, compute_g}, energy::compute_electronic_energy_expectation_value};



pub fn scf_cycle(s_mat: &Array2<f64>, t_mat: &Array2<f64>, vne: &Array2<f64>, vee: &Array4<f64>, tol: f64, max_iter: u8, mol: &Molecule) -> f64 {
    let mut electronic_energy = 0.0;

    let nbasis = mol.atoms.len();
    let mut dens_mat = Array2::<f64>::zeros((nbasis, nbasis));

    let (eigvals, eigvecs) = s_mat.clone().eigh(UPLO::Lower).unwrap();

    let mut lambda_inv_sqrt = Array2::<f64>::zeros((nbasis, nbasis));
    for i in 0..nbasis {
        lambda_inv_sqrt[(i, i)] = 1.0 / eigvals[i].sqrt();
    }

    let s_inv_sqrt = eigvecs.dot(&lambda_inv_sqrt.dot(&eigvecs.t()));

    for _step in 0..max_iter {
        let e_old = electronic_energy;

        let g = compute_g::compute_g(&dens_mat, vee);

        let f = t_mat + vne + &g;

        let f_units = s_inv_sqrt.dot(&f.dot(&s_inv_sqrt));

        let (_e_vals, e_vecs) = f_units.eigh(UPLO::Lower).unwrap();

        let mos = s_inv_sqrt.dot(&e_vecs);

        dens_mat = compute_density_matrix::compute_density_matrix(mos);

        electronic_energy = compute_electronic_energy_expectation_value::compute_electronic_energy_expectation_value(&dens_mat, &t_mat, &vne, &g);

        if (electronic_energy - e_old).abs() < tol {
            break;
        }
    }
    return electronic_energy;
}
