pub trait Potential {
    /// Derivada dU/dl: positiva en tensión, negativa en compresión.
    /// La fuerza sobre i apunta hacia j cuando este valor es positivo.
    fn force(&self, length: f64, initial_length: f64) -> f64;
    fn energy(&self, length: f64, initial_length: f64) -> f64;
}
