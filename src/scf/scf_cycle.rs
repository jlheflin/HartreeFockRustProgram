use ndarray::{Array1, Array2, Array4};
use ndarray_linalg::{Eigh, Solve, UPLO};

use crate::{classes::molecule::Molecule, density::{compute_density_matrix, compute_g}, energy::compute_electronic_energy_expectation_value};


struct DIISManager {
    fock_history: Vec<Array2<f64>>,
    error_history: Vec<Array2<f64>>,
    max_diis: usize,
}

impl DIISManager {
    pub fn new(max_diis: usize) -> Self {
        DIISManager { fock_history: Vec::new(), error_history: Vec::new(), max_diis }
    }

    pub fn add(&mut self, fock: Array2<f64>, error: Array2<f64>) {
        self.fock_history.push(fock);
        self.error_history.push(error);

        if self.fock_history.len() > self.max_diis {
            self.fock_history.remove(0);
            self.error_history.remove(0);
        }
    }

    pub fn extrapolate(&self) -> Option<Array2<f64>> {
        let n = self.error_history.len();
        if n == 0 {
            return None;
        }

        let mut b = Array2::<f64>::zeros((n + 1, n + 1));
        let mut rhs = Array1::<f64>::zeros(n + 1);
        rhs[n] = -1.0;

        for i in 0..n {
            for j in 0..n {
                let dot = (&self.error_history[i] * &self.error_history[j]).sum();
                b[(i, j)] = dot;
            }
            b[(i, n)] = -1.0;
            b[(n, i)] = -1.0;
        }

        let coeffs = match b.solve_into(rhs) {
            Ok(c) => c,
            Err(_) => return None,
        };

        let shape = self.fock_history[0].shape();

        let mut fock_new = Array2::<f64>::zeros((shape[0], shape[1]));
        for i in 0..n {
            fock_new.scaled_add(coeffs[i], &self.fock_history[i]);
        }

        Some(fock_new)
    }

    pub fn ready(&self) -> bool {
        self.fock_history.len() >= 2
    }
}


pub fn scf_cycle(s_mat: &Array2<f64>, t_mat: &Array2<f64>, vne: &Array2<f64>, vee: &Array4<f64>, tol: f64, max_iter: u8, mol: &Molecule, charge: i8, diis_enabled: bool) -> f64 {
    let mut electronic_energy = 0.0;

    let nbasis = mol.atoms.len();
    let mut dens_mat = Array2::<f64>::zeros((nbasis, nbasis));

    let (eigvals, eigvecs) = s_mat.clone().eigh(UPLO::Lower).unwrap();

    let mut lambda_inv_sqrt = Array2::<f64>::zeros((nbasis, nbasis));
    for i in 0..nbasis {
        lambda_inv_sqrt[(i, i)] = 1.0 / eigvals[i].sqrt();
    }

    let s_inv_sqrt = eigvecs.dot(&lambda_inv_sqrt.dot(&eigvecs.t()));

    let mut electrons: i8 = 0;
    for atom in mol.atoms.iter() {
        electrons += atom.atomic_number as i8;
    }

    electrons -= charge;
    let n_occ = electrons / 2;

    let mut diis = DIISManager::new(6);

    for _step in 0..max_iter {
        let e_old = electronic_energy;
        let g = compute_g::compute_g(&dens_mat, vee);
        let mut f = t_mat + vne + &g;

        if _step > 0 {
            let err = f.dot(&dens_mat).dot(s_mat) - s_mat.dot(&dens_mat).dot(&f);
            diis.add(f.clone(), err);

            if diis.ready() && diis_enabled {
                if let Some(f_extrapolated) = diis.extrapolate() {
                    f = f_extrapolated;
                }
            }
        }

        let f_units = s_inv_sqrt.dot(&f.dot(&s_inv_sqrt));

        let (_e_vals, e_vecs) = f_units.eigh(UPLO::Lower).unwrap();

        let mos = s_inv_sqrt.dot(&e_vecs);

        dens_mat = compute_density_matrix::compute_density_matrix(mos, n_occ);

        electronic_energy = compute_electronic_energy_expectation_value::compute_electronic_energy_expectation_value(&dens_mat, &t_mat, &vne, &g);

        if (electronic_energy - e_old).abs() < tol {
            break;
        }
    }
    return electronic_energy;
}
