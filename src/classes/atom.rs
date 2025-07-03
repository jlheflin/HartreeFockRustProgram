use nalgebra::Vector3;

use super::primitive_gaussian::PrimitiveGaussian;

pub struct Atom {
    pub atomic_number: u8,
    pub pgs: Vec<PrimitiveGaussian>
}

impl Atom {
    pub fn new(atomic_number: u8, pgs: Vec<PrimitiveGaussian>) -> Self {
        Self { atomic_number, pgs }
    }

    pub fn coords(&self) -> Vector3<f64> {
        self.pgs[0].coordinates
    }
}
