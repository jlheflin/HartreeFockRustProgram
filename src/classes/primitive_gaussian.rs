use nalgebra::Vector3;

pub struct PrimitiveGaussian {
    pub alpha: f64,
    pub coeff: f64,
    pub coordinates: Vector3<f64>
}

impl PrimitiveGaussian {
    pub fn new(alpha: f64, coeff: f64, coordinates: Vector3<f64>) -> Self {
        Self { alpha, coeff, coordinates }
    }

    pub fn a(&self) -> f64 {
        return (2.0 * self.alpha / std::f64::consts::PI).powf(0.75)
    }
}
