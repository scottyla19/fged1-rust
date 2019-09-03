extern crate game_math;

use game_math::vector::Vector3D;

fn main() {
    let a = Vector3D{x:1.0, y:2.0,z:3.0};
    let b = Vector3D{x:2.0, y:-1.0,z:5.0};
    println!("{:?}", a * b);
   
}