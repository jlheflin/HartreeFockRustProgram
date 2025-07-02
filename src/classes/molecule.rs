use nalgebra::{coordinates, Vector3};

use crate::classes::atom::Atom;

pub struct Molecule {
    pub atoms: Vec<Atom>,
}

impl Molecule {
    pub fn new() -> Self {
        Molecule { atoms: Vec::new() }
    }

    pub fn add_atom(&mut self, atom: Atom) {
        self.atoms.push(atom);
    }

    pub fn z_list(&self) -> Vec<u8> {
        self.atoms.iter().map(|atom| atom.atomic_number).collect()
    }

    pub fn coords_list(&self) -> Vec<Vector3<f64>> {
        self.atoms.iter().map(|atom| atom.coords()).collect()
    }
}
