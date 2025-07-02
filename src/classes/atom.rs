use nalgebra::Vector3;

use super::primitive_gaussian::PrimitiveGaussian;

pub struct Atom {
    pub atomic_number: u8,
    pub name: String,
    pub pgs: Vec<PrimitiveGaussian>
}

impl Atom {
    pub fn new(atomic_number: u8, name: String, pgs: Vec<PrimitiveGaussian>) -> Self {
        Self { atomic_number, name, pgs }
    }

    pub fn coords(&self) -> Vector3<f64> {
        self.pgs[0].coordinates
    }
}
