use std::{
    arch::x86_64::{
        _mm_add_ps, _mm_cvtss_f32, _mm_hadd_ps, _mm_insert_ps, _mm_mul_ps, _mm_set1_ps,
    },
    hint::unreachable_unchecked,
};

pub type Vec4 = core::arch::x86_64::__m128;
pub type Mat4 = [Vec4; 4];

/// Initialize a vector.
#[inline(always)]
pub const fn splat(x0: f32, x1: f32, x2: f32, x3: f32) -> Vec4 {
    #[repr(align(16))]
    struct Vec([f32; 4]);
    unsafe { std::mem::transmute(Vec([x0, x1, x2, x3]).0) }
}

/// Multiply vector by constant.
#[inline(always)]
pub fn cmul(c: f32, v: Vec4) -> Vec4 {
    unsafe {
        let c = _mm_set1_ps(c);
        _mm_mul_ps(c, v)
    }
}

/// Add two vectors.
#[inline(always)]
pub fn add(a: Vec4, b: Vec4) -> Vec4 {
    unsafe { _mm_add_ps(a, b) }
}

/// Store a variable at an index in a vector.
pub fn store1(v: Vec4, c: f32, n: u8) -> Vec4 {
    debug_assert!(n < 4);
    unsafe {
        let b = _mm_set1_ps(c);
        match n {
            0 => _mm_insert_ps::<0x00>(v, b),
            1 => _mm_insert_ps::<0x10>(v, b),
            2 => _mm_insert_ps::<0x20>(v, b),
            3 => _mm_insert_ps::<0x40>(v, b),
            _ => unreachable_unchecked(),
        }
    }
}

/// Dot product of two vectors.
#[inline(always)]
pub fn dot(a: Vec4, b: Vec4) -> f32 {
    unsafe {
        let mul = _mm_mul_ps(a, b);
        let add = _mm_hadd_ps(mul, mul);
        let dot = _mm_hadd_ps(add, add);
        _mm_cvtss_f32(dot)
    }
}
