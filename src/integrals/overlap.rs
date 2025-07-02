use crate::classes::molecule::Molecule;
use nalgebra::DMatrix;
use std::f64::consts::PI;

pub fn overlap(mol: &Molecule) -> DMatrix<f64> {
    let nbasis = mol.atoms.len();
    let mut s_mat = DMatrix::<f64>::zeros(nbasis, nbasis);
    let atoms = &mol.atoms;
    let coords = &mol.coords_list();

    for i in 0..nbasis {
        for j in 0..nbasis {
            let nprimitives_i = atoms[i].pgs.len();
            let prim_i = &atoms[i].pgs;
            let nprimitives_j = atoms[j].pgs.len();
            let prim_j = &atoms[j].pgs;

            for k in 0..nprimitives_i {
                for l in 0..nprimitives_j {
                    let n = prim_i[k].a() * prim_j[l].a();
                    let p = prim_i[k].alpha + prim_j[l].alpha;
                    let q = prim_i[k].alpha * prim_j[l].alpha / p;
                    let q_vec = coords[i] - coords[j];
                    let q2 = q_vec.dot(&q_vec);

                    s_mat[(i, j)] += n * prim_i[k].coeff * prim_j[l].coeff * (-q * q2).exp() * (PI / p).powf(1.5);
                    
                }
            }
        }
    }
    return s_mat;
}
