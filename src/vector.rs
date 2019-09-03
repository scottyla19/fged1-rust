
use std::ops::{Add, Sub, Mul,Div};

#[derive(Debug, Copy, Clone)]
pub struct Vector3D{
    pub x: f32,
    pub y: f32,
    pub z: f32,
}

impl Add<&Vector3D> for &Vector3D {
    type Output = Vector3D;
    fn add(self, v2: &Vector3D) -> Vector3D{
        Vector3D{x: self.x + v2.x, y:self.y + v2.y, z:self.z + v2.z}
    }
}

impl Sub<&Vector3D> for &Vector3D {
    type Output = Vector3D;
    fn sub(self, v2: &Vector3D) -> Vector3D{
        Vector3D{x: self.x - v2.x, y:self.y - v2.y, z:self.z - v2.z}
    }
}
// own and consume both vectors
impl Mul<Vector3D> for Vector3D {
    type Output = f32;
    fn mul(self, v2: Vector3D) -> f32{
        self.x * v2.x + self.y * v2.y + self.z * v2.z
    }
}
// borrow both vectors
impl Mul<&Vector3D> for &Vector3D {
    type Output = f32;
    fn mul(self, v2: &Vector3D) -> f32{
        self.x * v2.x + self.y * v2.y + self.z * v2.z
    }
}
// own and consume
impl Mul<f32> for Vector3D {
    type Output = Vector3D;
    fn mul(self, s: f32) -> Vector3D{
        Vector3D{x: self.x * s, y:self.y * s, z:self.z * s}
    }
}

impl Mul<f32> for &Vector3D {
    type Output = Vector3D;
    fn mul(self, s: f32) -> Vector3D{
        Vector3D{x: self.x * s, y:self.y * s, z:self.z * s}
    }
}
// 
impl Mul<&f32> for Vector3D {
    type Output = Vector3D;
    fn mul(self, s: &f32) -> Vector3D{
        Vector3D{x: self.x * s, y:self.y * s, z:self.z * s}
    }
}
impl Div<&f32> for Vector3D {
    type Output = Vector3D;
    fn div(self, s: &f32) -> Vector3D{
        Vector3D{x: self.x / s, y:self.y /s, z:self.z /s}
    }
}
impl Vector3D {
    pub fn dot<'a>(v1: &'a Vector3D, v2: &Vector3D) -> f32 {
        v1 * v2
    }
    pub fn cross<'a>(v1: &'a Vector3D, v2: &Vector3D) -> Vector3D {
        Vector3D {x: v1.y*v2.z - v1.z*v2.y, y: v1.z*v2.x - v1.x*v2.z, z: v1.x*v2.y - v1.y*v2.x}
    }
    pub fn project<'a>(v1: &'a Vector3D, v2: &Vector3D) -> Vector3D {
        v2*((v1*v2)/(v2*v2))
    }

    pub fn reject<'a>(v1: &'a Vector3D, v2: &Vector3D) -> Vector3D {
        let proj = Vector3D::project(v1, v2);
        Vector3D{x:v1.x-proj.x, y:v1.y - proj.y, z:v1.z - proj.z}
    }

    pub fn normalize(&self) -> Vector3D {
        let mag = self.magnitude();
        Vector3D{x:self.x/mag, y:self.y/mag, z:self.z/mag}
    }
    pub fn magnitude(&self) -> f32 {
        (self.x*self.x + self.y*self.y + self.z*self.z).sqrt()
    }
}


#[cfg(test)]

mod tests {
    use crate::vector::Vector3D;
    const V1 : Vector3D = Vector3D{x:2.0,y:3.0,z:4.0};
    const V2 : Vector3D = Vector3D{x:5.0,y:6.0,z:7.0};
    const V3 : Vector3D = Vector3D{x:3.0,y:4.0,z:0.0};
    #[test]
    fn test_add_two_vector3ds() {
        let sum = &V1 + &V2;
        assert_eq!(sum.x, 7.0);
        assert_eq!(sum.y, 9.0);
        assert_eq!(sum.z, 11.0);
    }
    #[test]
    fn test_sub_two_vector3ds() {
        let subt = &V1 - &V2;
        assert_eq!(subt.x, -3.0);
        assert_eq!(subt.y, -3.0);
        assert_eq!(subt.z, -3.0);
    }
    #[test]
    fn test_mul_two_vector3ds() {
        let dot = &V1 * &V2;
        assert_eq!(dot, 56.0);
    }
    
     #[test]
    fn test_mul_scalar() {
        let prod = V1 * 2.0;
        assert_eq!(prod.x, 4.0);
        assert_eq!(prod.y, 6.0);
        assert_eq!(prod.z, 8.0);
    }
    #[test]
    fn test_borrow_mul_scalar() {
        let scalar = 2.0;
        let prod = V1 * &scalar;
        assert_eq!(prod.x, 4.0);
        assert_eq!(prod.y, 6.0);
        assert_eq!(prod.z, 8.0);
    }
     #[test]
    fn test_div_scalar() {
        let scalar = 2.0;
        let quot = V1 / &scalar;
        assert_eq!(quot.x, 1.0);
        assert_eq!(quot.y, 1.5);
        assert_eq!(quot.z, 2.0);
    }
    #[test]
    fn test_dot_two_vector3ds() {
        let dot = Vector3D::dot(&V1, &V2);
        assert_eq!(dot, 56.0);
    }
    #[test]
    fn test_project_two_vector3ds() {
        let proj = Vector3D::project(&V1, &V2);
        assert_eq!(proj.x, 28.0/11.0);
        assert_eq!(proj.y, 168.0/55.0);
        assert_eq!(proj.z, 196.0/55.0);
    }

    #[test]
    fn test_reject_two_vector3ds() {
        let rej = Vector3D::reject(&V1, &V2);
        assert_eq!(rej.x, 2.0 - (28.0/11.0));
        assert_eq!(rej.y, 3.0 -(168.0/55.0));
        assert_eq!(rej.z, 4.0 - (196.0/55.0));
    }
    #[test]
    fn test_cross_two_vector3ds() {
        let cross = Vector3D::cross(&V1, &V2);
        assert_eq!(cross.x, -3.0);
        assert_eq!(cross.y, 6.0);
        assert_eq!(cross.z, -3.0);
    }
    #[test]
    fn test_magnitude(){
        let mag = V1.magnitude();
        assert_eq!(29.0_f32.sqrt(), mag);
    }
    #[test]
    fn test_normalize(){
        let norm = V1.normalize();
        let mag = V1.magnitude();
        assert_eq!(V1.x/mag, norm.x);
        assert_eq!(V1.y/mag, norm.y);
        assert_eq!(V1.z/mag, norm.z);
    }
}
