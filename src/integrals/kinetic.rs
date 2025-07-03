use crate::classes::molecule::Molecule;
use ndarray::Array2;
use std::f64::consts::PI;

pub fn kinetic(mol: &Molecule) -> Array2<f64> {
    let nbasis = mol.atoms.len();
    let mut t_mat = Array2::<f64>::zeros((nbasis, nbasis));
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
                    let cacb = prim_i[k].coeff * prim_j[l].coeff;
                    let p = prim_i[k].alpha + prim_j[l].alpha;
                    let p_coords = prim_i[k].alpha * coords[i] + prim_j[l].alpha * coords[j];
                    let pp = p_coords / p;
                    let pg = pp - coords[j];
                    let pgx2 = pg[0] * pg[0];
                    let pgy2 = pg[1] * pg[1];
                    let pgz2 = pg[2] * pg[2];

                    let q = prim_i[k].alpha * prim_j[l].alpha / p;
                    let q_vec = coords[i] - coords[j];
                    let q2 = q_vec.dot(&q_vec);

                    let s = (-q * q2).exp() * (PI / p).powf(1.5) * n * cacb;

                    t_mat[(i, j)] += 3.0 * prim_j[l].alpha * s;
                    t_mat[(i, j)] -= 2.0 * prim_j[l].alpha * prim_j[l].alpha * s * (pgx2 + 0.5 / p);
                    t_mat[(i, j)] -= 2.0 * prim_j[l].alpha * prim_j[l].alpha * s * (pgy2 + 0.5 / p);
                    t_mat[(i, j)] -= 2.0 * prim_j[l].alpha * prim_j[l].alpha * s * (pgz2 + 0.5 / p);
                    
                }
            }
        }
    }
    return t_mat;
}
