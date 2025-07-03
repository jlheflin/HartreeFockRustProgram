use std::f64::consts::PI;
use ndarray::Array4;
use crate::{classes::molecule::Molecule, integrals::boys::{self}};



pub fn electron_electron_repulsion(mol: &Molecule) -> Array4<f64> {
    let nbasis = mol.atoms.len();
    let atoms = &mol.atoms;

    let mut v_ee = Array4::<f64>::zeros((nbasis, nbasis, nbasis, nbasis));

    for i in 0..nbasis {
        for j in 0..nbasis {
            for k in 0..nbasis {
                for l in 0..nbasis {

                    let nprimitives_i = atoms[i].pgs.len();
                    let prim_i = &atoms[i].pgs;
                    let nprimitives_j = atoms[j].pgs.len();
                    let prim_j = &atoms[j].pgs;
                    let nprimitives_k = atoms[k].pgs.len();
                    let prim_k = &atoms[k].pgs;
                    let nprimitives_l = atoms[l].pgs.len();
                    let prim_l = &atoms[l].pgs;

                    for ii in 0..nprimitives_i {
                        for jj in 0..nprimitives_j {
                            for kk in 0..nprimitives_k {
                                for ll in 0..nprimitives_l {

                                    let n = prim_i[ii].a() * prim_j[jj].a() * prim_k[kk].a() * prim_l[ll].a();
                                    let c = prim_i[ii].coeff * prim_j[jj].coeff * prim_k[kk].coeff * prim_l[ll].coeff;
                                    let pij = prim_i[ii].alpha + prim_j[jj].alpha;
                                    let pkl = prim_k[kk].alpha + prim_l[ll].alpha;

                                    let pij_coords = prim_i[ii].alpha * prim_i[ii].coordinates + prim_j[jj].alpha * prim_j[jj].coordinates;
                                    let pkl_coords = prim_k[kk].alpha * prim_k[kk].coordinates + prim_l[ll].alpha * prim_l[ll].coordinates;

                                    let ppij = pij_coords / pij;
                                    let ppkl = pkl_coords / pkl;

                                    let ppijppkl = ppij - ppkl;
                                    let ppijppkl2 = ppijppkl.dot(&ppijppkl);

                                    let denom = 1.0 / pij + 1.0 / pkl;

                                    let qij = prim_i[ii].alpha * prim_j[jj].alpha / pij;
                                    let qkl = prim_k[kk].alpha * prim_l[ll].alpha / pkl;

                                    let qij_coords = prim_i[ii].coordinates - prim_j[jj].coordinates;
                                    let qkl_coords = prim_k[kk].coordinates - prim_l[ll].coordinates;

                                    let q2ij = qij_coords.dot(&qij_coords);
                                    let q2kl = qkl_coords.dot(&qkl_coords);

                                    let term1 = 2.0 * PI.powi(2) / (pij * pkl);
                                    let term2 = (PI / (pij + pkl)).sqrt();
                                    let term3 = (-qij * q2ij).exp();
                                    let term4 = (-qkl * q2kl).exp();
                                    let b = boys::boys(ppijppkl2/denom, 0.0);

                                    v_ee[(i, j, k, l)] += n * c * term1 * term2 * term3 * term4 * b;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return v_ee;
}
