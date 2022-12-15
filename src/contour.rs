use std::collections::{HashMap, HashSet};

use crate::{Triangulation, Point2, HasPosition};
trait HasValue
{
    fn get_value() -> f64;
}

// fn contour<T>(mesh: &T, levels: Vec<f64>) -> Vec<Point2<<<T as Triangulation>::Vertex as HasPosition>::Scalar>> where T: Triangulation + {<T as Triangulation>::Vertex: HasValue}: HasValue
// {
//     let map: HashMap<f64, T::Face> = HashMap::new();
//     let mut levels_copy: Vec<f64> = levels.into_iter().collect();
//     for face in mesh.inner_faces()
//     {
//         let v1 = face.vertices()[0].data();
//         let v2 = face.vertices()[0].data();
//         let v3 = face.vertices()[0].data();
        
//     }
//     return Vec::new();
// }
    
    
