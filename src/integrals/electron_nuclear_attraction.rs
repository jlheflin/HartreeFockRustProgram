use std::f64::consts::PI;
use nalgebra::{Vector3};
use ndarray::Array2;

use crate::{classes::molecule::Molecule, integrals::boys};



pub fn electron_nuclear_attraction(mol: &Molecule) -> Array2<f64> {
    let natoms = mol.atoms.len();
    let nbasis = mol.atoms.len();
    let atoms = &mol.atoms;
    let coords = &mol.coords_list();
    let z_list = mol.z_list();

    let mut v_ne = Array2::<f64>::zeros((nbasis, nbasis));

    let mut unique_coords: Vec<Vector3<f64>> = Vec::new();
    let eps = 1e-8;

    for atom in mol.atoms.iter() {
        for pg in atom.pgs.iter() {
            let coord = pg.coordinates;

            if !unique_coords.iter().any(|c| (c - coord).norm() < eps) {
                unique_coords.push(coord);
            }
        }
    }

    for atom in 0..natoms {
        for i in 0..nbasis {
            for j in 0..nbasis {
                let nprimitives_i = atoms[i].pgs.len();
                let prim_i = &atoms[i].pgs;
                let nprimitives_j = atoms[j].pgs.len();
                let prim_j = &atoms[j].pgs;

                for k in 0..nprimitives_i {
                    for l in 0..nprimitives_j {
                        
                        let n = prim_i[k].a() * prim_j[l].a();
                        let cacb = prim_i[k].coeff * prim_j[l].coeff;
                        let p = prim_i[k].alpha + prim_j[l].alpha;
                        let p_coords = prim_i[k].alpha * coords[i] + prim_j[l].alpha * coords[j];
                        let pp = p_coords / p;
                        let pg = pp - coords[atom];

                        let pg2 = pg.dot(&pg);

                        let q = prim_i[k].alpha * prim_j[l].alpha / p;
                        let q_vec = coords[i] - coords[j];
                        let q2 = q_vec.dot(&q_vec);

                        v_ne[(i, j)] += n * cacb * (-1.0 * z_list[atom] as f64) * (2.0 * PI / p) * (-q * q2).exp() * boys::boys(p * pg2, 0.);
                    }
                }
            }
        }
    }
    return v_ne;
}
