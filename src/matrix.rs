
use crate::vector::Vector3D;
use std::fmt;
use std::ops::{Add, Mul, Sub};

// h/t https://stackoverflow.com/questions/49037111/alternatives-to-matching-floating-point-ranges
trait InRange {
    fn in_range(&self, begin: Self, end: Self) -> bool;
}

impl InRange for f32 {
    fn in_range(&self, begin: f32, end: f32) -> bool {
        *self >= begin && *self < end
    }
}

#[derive(Debug, Copy, Clone)]
pub struct Matrix3D {
    pub r1: Vector3D,
    pub r2: Vector3D,
    pub r3: Vector3D,
}
impl Add<&Matrix3D> for &Matrix3D {
    type Output = Matrix3D;
    fn add(self, m2: &Matrix3D) -> Matrix3D {
        Matrix3D {
            r1: &self.r1 + &m2.r1,
            r2: &self.r2 + &m2.r2,
            r3: &self.r3 + &m2.r3,
        }
    }
}
impl Sub<&Matrix3D> for &Matrix3D {
    type Output = Matrix3D;
    fn sub(self, m2: &Matrix3D) -> Matrix3D {
        Matrix3D {
            r1: &self.r1 - &m2.r1,
            r2: &self.r2 - &m2.r2,
            r3: &self.r3 - &m2.r3,
        }
    }
}
impl Mul<&Matrix3D> for &Matrix3D {
    type Output = Matrix3D;
    fn mul(self, m2: &Matrix3D) -> Matrix3D {
        let b_trans = m2.transpose();
        let new_r1 = Vector3D {
            x: &self.r1 * &b_trans.r1,
            y: &self.r1 * &b_trans.r2,
            z: &self.r1 * &b_trans.r3,
        };
        let new_r2 = Vector3D {
            x: &self.r2 * &b_trans.r1,
            y: &self.r2 * &b_trans.r2,
            z: &self.r2 * &b_trans.r3,
        };
        let new_r3 = Vector3D {
            x: &self.r3 * &b_trans.r1,
            y: &self.r3 * &b_trans.r2,
            z: &self.r3 * &b_trans.r3,
        };
        Matrix3D {
            r1: new_r1,
            r2: new_r2,
            r3: new_r3,
        }
    }
}
impl Mul<&Vector3D> for &Matrix3D {
    type Output = Vector3D;
    fn mul(self, v: &Vector3D) -> Vector3D {
        Vector3D {
            x: &self.r1 * v,
            y: &self.r2 * v,
            z: &self.r3 * v,
        }
    }
}
impl Mul<&f32> for &Matrix3D {
    type Output = Matrix3D;
    fn mul(self, s: &f32) -> Matrix3D {
        Matrix3D {
            r1: self.r1 * s,
            r2: self.r2 * s,
            r3: self.r3 * s,
        }
    }
}
impl Matrix3D {
    pub fn new_constant(value: f32) -> Matrix3D {
        Matrix3D {
            r1: Vector3D {
                x: value,
                y: value,
                z: value,
            },
            r2: Vector3D {
                x: value,
                y: value,
                z: value,
            },
            r3: Vector3D {
                x: value,
                y: value,
                z: value,
            },
        }
    }
    pub fn new_identity() -> Matrix3D {
        Matrix3D {
            r1: Vector3D {
                x: 1.0,
                y: 0.0,
                z: 0.0,
            },
            r2: Vector3D {
                x: 0.0,
                y: 1.0,
                z: 0.0,
            },
            r3: Vector3D {
                x: 0.0,
                y: 0.0,
                z: 1.0,
            },
        }
    }
    pub fn transpose(&self) -> Matrix3D {
        let new_r1 = Vector3D {
            x: self.r1.x,
            y: self.r2.x,
            z: self.r3.x,
        };
        let new_r2 = Vector3D {
            x: self.r1.y,
            y: self.r2.y,
            z: self.r3.y,
        };
        let new_r3 = Vector3D {
            x: self.r1.z,
            y: self.r2.z,
            z: self.r3.z,
        };
        Matrix3D {
            r1: new_r1,
            r2: new_r2,
            r3: new_r3,
        }
    }
    pub fn determinant(&self) -> f32 {
        let col_vectors = self.transpose();
        let a = col_vectors.r1;
        let b = col_vectors.r2;
        let c = col_vectors.r3;
        a.x * b.y * c.z + a.y * b.z * c.x + a.z * b.x * c.y
            - a.x * b.z * c.y
            - a.y * b.x * c.z
            - a.z * b.y * c.x
    }
    pub fn inverse(&self) -> Result<Matrix3D, &str> {
        let det = self.determinant();
        // let det = det as i32;
        match det {
            x if x.in_range(-0.00001, 0.00001) => {
                Err("Determinant = 0. Cannot calculate the inverse.")
            }
            _ => Ok(self.calc_inverse()),

        }
    }
    fn calc_inverse(&self) -> Matrix3D {
        let cols = self.transpose();
        let inv_det = 1.0 / self.determinant();
        let a = cols.r1;
        let b = cols.r2;
        let c = cols.r3;
        let new_r1 = Vector3D::cross(&b, &c);
        let new_r2 = Vector3D::cross(&c, &a);
        let new_r3 = Vector3D::cross(&a, &b);
        Matrix3D {
            r1: new_r1 * inv_det,
            r2: new_r2 * inv_det,
            r3: new_r3 * inv_det,
        }
    }

}

impl fmt::Display for Matrix3D {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "\n{}  {}  {}\n{}  {}  {}\n{}  {}  {}\n",
            self.r1.x,
            self.r1.y,
            self.r1.z,
            self.r2.x,
            self.r2.y,
            self.r2.z,
            self.r3.x,
            self.r3.y,
            self.r3.z,
        )
    }
}


#[cfg(test)]

mod tests {
    use crate::matrix::Matrix3D;
    use crate::vector::Vector3D;
    const M0: Matrix3D = Matrix3D {
        r1: Vector3D {
            x: -9.0,
            y: -8.0,
            z: -7.0,
        },
        r2: Vector3D {
            x: -6.0,
            y: -5.0,
            z: -4.0,
        },
        r3: Vector3D {
            x: -3.0,
            y: -2.0,
            z: -1.0,
        },
    };
    const M1: Matrix3D = Matrix3D {
        r1: Vector3D {
            x: 1.0,
            y: 2.0,
            z: 3.0,
        },
        r2: Vector3D {
            x: 4.0,
            y: 5.0,
            z: 6.0,
        },
        r3: Vector3D {
            x: 7.0,
            y: 8.0,
            z: 9.0,
        },
    };
    const M2: Matrix3D = Matrix3D {
        r1: Vector3D {
            x: 1.0,
            y: 2.0,
            z: 4.0,
        },
        r2: Vector3D {
            x: 4.0,
            y: 5.0,
            z: 6.0,
        },
        r3: Vector3D {
            x: 7.0,
            y: 8.0,
            z: 9.0,
        },
    };
    #[test]
    fn test_new_constant_matrix() {
        let constant = 0.0;
        let m = Matrix3D::new_constant(constant);
        assert_eq!(m.r1.x, constant);
        assert_eq!(m.r1.y, constant);
        assert_eq!(m.r1.z, constant);
        assert_eq!(m.r2.x, constant);
        assert_eq!(m.r2.y, constant);
        assert_eq!(m.r2.z, constant);
        assert_eq!(m.r3.x, constant);
        assert_eq!(m.r3.y, constant);
        assert_eq!(m.r3.z, constant);

    }
    #[test]
    fn test_new_identity_matrix() {
        let m = Matrix3D::new_identity();
        assert_eq!(m.r1.x, 1.0);
        assert_eq!(m.r1.y, 0.0);
        assert_eq!(m.r1.z, 0.0);
        assert_eq!(m.r2.x, 0.0);
        assert_eq!(m.r2.y, 1.0);
        assert_eq!(m.r2.z, 0.0);
        assert_eq!(m.r3.x, 0.0);
        assert_eq!(m.r3.y, 0.0);
        assert_eq!(m.r3.z, 1.0);
    }
    #[test]
    fn test_add_matrices() {
        let m = &M0 + &M1;
        assert_eq!(m.r1.x, -8.0);
        assert_eq!(m.r1.y, -6.0);
        assert_eq!(m.r1.z, -4.0);
        assert_eq!(m.r2.x, -2.0);
        assert_eq!(m.r2.y, 0.0);
        assert_eq!(m.r2.z, 2.0);
        assert_eq!(m.r3.x, 4.0);
        assert_eq!(m.r3.y, 6.0);
        assert_eq!(m.r3.z, 8.0);
    }
    #[test]
    fn test_sub_matrices() {
        let m = &M0 - &M1;
        assert_eq!(m.r1.x, -10.0);
        assert_eq!(m.r1.y, -10.0);
        assert_eq!(m.r1.z, -10.0);
        assert_eq!(m.r2.x, -10.0);
        assert_eq!(m.r2.y, -10.0);
        assert_eq!(m.r2.z, -10.0);
        assert_eq!(m.r3.x, -10.0);
        assert_eq!(m.r3.y, -10.0);
        assert_eq!(m.r3.z, -10.0);
    }
    #[test]
    fn test_transpose_matrix() {
        let m = &M1.transpose();
        assert_eq!(m.r1.x, 1.0);
        assert_eq!(m.r1.y, 4.0);
        assert_eq!(m.r1.z, 7.0);
        assert_eq!(m.r2.x, 2.0);
        assert_eq!(m.r2.y, 5.0);
        assert_eq!(m.r2.z, 8.0);
        assert_eq!(m.r3.x, 3.0);
        assert_eq!(m.r3.y, 6.0);
        assert_eq!(m.r3.z, 9.0);
    }
    #[test]
    fn test_multiply_matrices() {
        let m = &M0 * &M1;
        assert_eq!(m.r1.x, -90.0);
        assert_eq!(m.r1.y, -114.0);
        assert_eq!(m.r1.z, -138.0);
        assert_eq!(m.r2.x, -54.0);
        assert_eq!(m.r2.y, -69.0);
        assert_eq!(m.r2.z, -84.0);
        assert_eq!(m.r3.x, -18.0);
        assert_eq!(m.r3.y, -24.0);
        assert_eq!(m.r3.z, -30.0);
    }
    #[test]
    fn test_scalar_multiply_matrices() {
        let s = 2.0;
        let m = &M1 * &s;
        assert_eq!(m.r1.x, 2.0);
        assert_eq!(m.r1.y, 4.0);
        assert_eq!(m.r1.z, 6.0);
        assert_eq!(m.r2.x, 8.0);
        assert_eq!(m.r2.y, 10.0);
        assert_eq!(m.r2.z, 12.0);
        assert_eq!(m.r3.x, 14.0);
        assert_eq!(m.r3.y, 16.0);
        assert_eq!(m.r3.z, 18.0);
    }
    #[test]
    fn test_vector_multiply_matrices() {
        let v = Vector3D {
            x: 2.0,
            y: 3.0,
            z: 4.0,
        };
        let prod = &M1 * &v;
        assert_eq!(prod.x, 20.0);
        assert_eq!(prod.y, 47.0);
        assert_eq!(prod.z, 74.0);
    }
    #[test]
    fn test_matrix_determinant() {
        let det = M2.determinant();
        assert_eq!(det, -3.0);
    }
    #[test]
    fn test_invertable_matrix() {
        let inv = M2.inverse().unwrap();
        println!("{}", inv);
        assert_eq!(inv.r1.x, 1.0);
        assert_eq!(format!("{:.5}", inv.r1.y), format!("{:.5}", -14.0 / 3.0));
        assert_eq!(inv.r1.z, 8.0 / 3.0);
        assert_eq!(inv.r2.x, -2.0);
        assert_eq!(inv.r2.y, 19.0 / 3.0);
        assert_eq!(format!("{:.5}", inv.r2.z), format!("{:.5}", -10.0 / 3.0));
        assert_eq!(inv.r3.x, 1.0);
        assert_eq!(inv.r3.y, -2.0);
        assert_eq!(inv.r3.z, 1.0);
    }
    #[test]
    fn test_non_invertable_matrix() {
        let inv = M1.inverse();
        assert!(inv.is_err());
        // assert_eq!(inv., "Determinant = 0. Cannot calculate the inverse." );
    }
}