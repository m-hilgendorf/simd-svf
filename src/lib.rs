pub mod arch;

pub struct SVF {
    a: [arch::Vec4; 3],
    b: [arch::Vec4; 2],
    s: arch::Vec4,
}

impl Default for SVF {
    fn default() -> Self {
        Self::new()
    }
}

impl SVF {
    /// Create a new filter.
    pub fn new() -> Self {
        Self {
            a: [arch::zero(), arch::zero(), arch::zero()],
            b: [arch::zero(), arch::zero()],
            s: arch::zero(),
        }
    }

    /// Set the internals to a particular frequency/Q.
    pub fn set(&mut self, normalized_frequency: f32, q_factor: f32) {
        // Validate inputs.
        debug_assert!((0.0..1.0).contains(&normalized_frequency));
        debug_assert!((0.5..10.0).contains(&q_factor));

        // Compute coefficients.
        let fb = q_factor.recip();
        let ff = (core::f32::consts::FRAC_PI_2 * normalized_frequency).tan();
        let gain = (1.0 + fb * ff + ff * ff).recip();

        // Compute matrices.
        let a = [
            arch::cmul(gain, arch::splat(1.0, -(ff + fb), -1.0, 0.0)),
            arch::cmul(gain, arch::splat(ff, 1.0, -ff, 0.0)),
            arch::cmul(gain, arch::splat(ff * ff, ff, 1.0 + fb, 0.0)),
        ];
        let b = [
            arch::add(arch::cmul(ff, a[0]), a[1]),
            arch::add(arch::cmul(ff, a[1]), a[2]),
        ];
        self.a = a;
        self.b = b;
    }

    /// Reset filter state.
    pub fn reset(&mut self) {
        self.s = arch::zero()
    }

    /// Compute highpass/bandpass/lowpass outputs for a given input.
    pub fn eval(&mut self, x: f32) -> [f32; 3] {
        // Concat input with the state vector
        let z = arch::store1(self.s, x, 0);

        // Compute outputs.
        let yhp = arch::dot(self.a[0], z);
        let ybp = arch::dot(self.a[1], z);
        let ylp = arch::dot(self.a[2], z);

        // Update state.
        let s1 = arch::dot(self.b[0], z);
        let s2 = arch::dot(self.b[1], z);
        self.s = arch::splat(0.0, s1, s2, 0.0);

        // Return output.
        [yhp, ylp, ybp]
    }
}
