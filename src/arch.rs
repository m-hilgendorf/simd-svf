#[cfg(target_arch = "x86_64")]
mod x86_64;

#[cfg(target_arch = "x86_64")]
pub use x86_64::*;

pub fn zero() -> Vec4 {
    splat(0.0, 0.0, 0.0, 0.0)
}
