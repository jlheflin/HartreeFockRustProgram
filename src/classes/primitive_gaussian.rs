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

    pub fn x(&self) -> f64 {
        self.coordinates[0]
    }

    pub fn y(&self) -> f64 {
        self.coordinates[1]
    }

    pub fn z(&self) -> f64 {
        self.coordinates[2]
    }

    pub fn a(&self) -> f64 {
        return (2.0 * self.alpha / std::f64::consts::PI).powf(0.75)
    }
}
