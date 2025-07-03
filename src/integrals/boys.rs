use rgsl::gamma_beta::{gamma, incomplete_gamma};

pub fn boys(x: f64, n: f64) -> f64 {
    if x.abs() < 1e-8 {
        return 1.0 / (2.0 * n + 1.0)
    }
    else {
        let a = n + 0.5;
        let gamma_p = incomplete_gamma::gamma_inc_P(a, x);
        return gamma_p * gamma::gamma(a) * (1.0 / (2.0 * x.powf(a)));
    }
}
