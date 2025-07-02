use crate::classes::molecule::Molecule;

pub fn nuclear_nuclear_repusion_energy(mol: &Molecule) -> f64 {

    let natoms = mol.z_list().len();
    let zlist = mol.z_list();
    let coords = mol.coords_list();

    let mut e_nn = 0.0;

    for i in 0..natoms {
        let zi = zlist[i];
        for j in (i + 1)..natoms {
            let zj = zlist[j];

            let rij = coords[i] - coords[j];
            let distance = rij.norm();


            e_nn += (zi as f64 * zj as f64) / distance;
        }
    }
    return e_nn;
}
