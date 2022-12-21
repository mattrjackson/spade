use std::collections::{HashMap, HashSet};
use std::convert::From;
use std::ops::Mul;
use num_traits::ToPrimitive;

use crate::handles::{DirectedEdgeHandle, UndirectedEdgeHandle, FixedDirectedEdgeHandle};
use crate::internals::DynamicHandleImpl;
use crate::{Triangulation, Point2, HasPosition, SpadeNum, DelaunayTriangulation};
use crate::delaunay_core::{FaceTag, InnerTag, bulk_load};

#[derive(Clone)]
pub struct ContourResult<S>
{
    pub level: f64,
    pub points: Vec<Point2<S>>
}
impl<S> ContourResult<S>
{
    pub fn new(level: f64) -> ContourResult<S>
    {
        ContourResult { level: level, points: Vec::new() }
    }
}

///
/// Defines the Interval Bounding State 
/// There are three states:
/// AtVertex: Contour level found at an edge vertex
/// Inside:   Contour level is on edge between the edge vertices
/// Outside:  Contour level outside the bounds for the edge
/// 
#[derive(PartialEq)]
enum IntervalState
{
    AtVertex = 0,
    Inside = 1,
    Outside = 2
}
///
/// Determines the interval state for a given value and bounding
/// values on either vertex of an edge.
/// 
fn in_interval(val:f64, val1:f64, val2: f64) -> IntervalState
{
    let r = (val1-val)*(val2-val);
    if r > 0.0
    {
        IntervalState::Outside
    }
    else if r < 0.0
    {
        IntervalState::Inside
    }
    else
    {
        IntervalState::AtVertex
    }
} 

///
/// Trait that pulls a single f64 from a given Vertex data structure.
/// 
pub trait GetValue
{
    fn get_value(&self) -> f64;
}
pub trait Contour : Triangulation where
    <Self as Triangulation>::Vertex : HasPosition + GetValue
{
    ///
    /// Processes a given face that has an edge that bounds a given contour level
    /// 
    fn process_face_for_edge(&self, val:f64, fixed_edge: &FixedDirectedEdgeHandle) ->
        Option<FixedDirectedEdgeHandle>
    {
        let mut result: Option<FixedDirectedEdgeHandle> = None;
        let edge = self.directed_edge(*fixed_edge);
        let face = edge.face().as_inner();
        if face.is_none() // we've reached an external face...
        {
            return None;
        }
        for e in edge.face().as_inner().unwrap().adjacent_edges()
        {
            if edge == e // don't return the input edge
            {
                continue;
            }     
            let edge_from = e.from();
            let edge_to = e.to();
            let v1 = edge_from.data();
            let v2 = edge_to.data();
            let estate = in_interval(val, v1.get_value(), v2.get_value());
            if estate == IntervalState::Inside
            {   
                if e.rev().is_outer_edge() // skip outer edges
                {
                    continue;
                }
                else
                {
                    return Some(e.fix());
                }
            }
            else if estate == IntervalState::AtVertex
            {
                // this value will only get returned if the edge opposite does not form the second point of the contour.
                result = Some(e.fix());

            }
        }
        return result;    
    }

    ///
    /// Given a starting edge, this routine will complete the computation of a contour.
    /// 
    fn complete_contour(&self, val:f64, edge: &FixedDirectedEdgeHandle, 
    results: &mut Vec<FixedDirectedEdgeHandle>, visited_faces: &mut HashSet<usize>)
    {
        let first_edge = edge.clone();
        let face = &self.directed_edge(*edge).face().as_inner().unwrap();    
        let mut edges: Vec<FixedDirectedEdgeHandle> = Vec::new();
        let mut current_edge = Some(first_edge);    

        if visited_faces.contains(&face.fix().index())
        {
            // This face has already been visited. Don't try to add it's edges again.
            return;
        }
        // tag this face as visited
        visited_faces.insert(face.fix().index());    
        // Wind through faces until we reach our original point or hit an outer boundary.
        loop
        {
            let new_edge = self.process_face_for_edge(val, &current_edge.unwrap());       
            current_edge = if new_edge.is_none() { None} else { Some(new_edge.unwrap())};
            if current_edge.is_some() 
            {
                // tag this edge's face as visited...
                let de = self.directed_edge(current_edge.unwrap());            
                visited_faces.insert(de.face().fix().index());
                // add edge to our list we will return...
                edges.push(current_edge.unwrap());            
                if current_edge.unwrap() == first_edge.rev()
                {
                    // if we've reached our original edge, the contour is closed.                
                    break;
                }
                // finally get the reversed edge and evaluate its face in the next iteration.
                current_edge = Some(current_edge.unwrap().rev());
            }
            else
            { 
                break;
            }
        }
        if current_edge.is_none() && !self.directed_edge(first_edge).rev().face().is_outer()
        {       
            // If we've arrived here, we have reached an outer edge but have potentially not finished drawing our contour.
            // So we reverse our search direction relative the first edge, compute the intersecting edges, and finally
            // return the completed contour.

            let mut edges2: Vec<FixedDirectedEdgeHandle> = Vec::new();
            // calling complete_contour with first_edge reversed.        
            self.complete_contour(val, &first_edge.rev(), &mut edges2, visited_faces);
            // since the previous calculation was in the opposite direction, we need to reverse the order of our edges...
            edges2.reverse();
            // now copy reordered edges to our edges vector...
            for e in &edges
            {
                edges2.push(e.clone());            
            }
            edges = edges2
        }    
        else if self.directed_edge(first_edge).rev().is_outer_edge() // if we got lucky and started at one boundary and ended at another, we should make sure to tag that last face so we don't try to iterate through it again...
        {
            visited_faces.insert(self.directed_edge(first_edge).rev().face().fix().index());
        }
        *results = edges;

    }

    fn compute_point(&self, level:f64, fixed_edge: &FixedDirectedEdgeHandle, take_log_data: bool) -> Point2<<Self::Vertex as HasPosition>::Scalar>
    {
        let edge = self.directed_edge(*fixed_edge);
        let x0 = edge.from().position().x;
        let y0 = edge.from().position().y;
        let x1 = edge.to().position().x;
        let y1 = edge.to().position().y;
        
        let mut v0:f64 = (*edge.from().data()).get_value();    
        let mut v1:f64 = (*edge.to().data()).get_value();
        let mut v = level.clone();
        if take_log_data
        {
            v0 = v0.ln();
            v1 = v1.ln();
            v = v.ln();
        }
        let f:<Self::Vertex as HasPosition>::Scalar = Into::into(((v-v0) / (v1-v0)) as f32);
        // println!("x1={:?} {}, x2={:?} {} ", edge.from().position(), edge.from().data().get_value(), edge.to().position(),edge.to().data().get_value());
        // println!("{:?} {:?}", Point2 { x:  x0 + (x1-x0)*f, y: y0 + (y1-y0)*f }, f);
        Point2 { x:  x0 + (x1-x0)*f, y: y0 + (y1-y0)*f }
    }

    fn contour(&self, level: f64, take_log_data: bool) -> Vec<ContourResult<<Self::Vertex as HasPosition>::Scalar>>
    {
        let mut visited_faces: HashSet<usize> = HashSet::new();
        let mut results = Vec::new();
        for face in self.inner_faces()
        {
            if visited_faces.contains(&face.fix().index())
            {
                continue;
            }
            for e in face.adjacent_edges()
            {                
                let val1:f64 = e.from().data().get_value();
                let val2:f64 = e.to().data().get_value();            
                if in_interval(level, val1, val2) == IntervalState::Inside
                {
                    
                    //println!("{:?} {:?} {:?} {} {} {}", face.vertices()[0].position(), face.vertices()[1].position(), face.vertices()[2].position(), face.vertices()[0].data().get_value(),face.vertices()[1].data().get_value(),face.vertices()[2].data().get_value());
                    let mut vec= Vec::new();                    
                    self.complete_contour(level, &e.fix(), &mut vec, &mut visited_faces);
                    let mut contour: ContourResult<<Self::Vertex as HasPosition>::Scalar> = ContourResult::new(level);
                    for e in vec
                    {          
                        let pt = self.compute_point(level, &e, take_log_data);
                        if contour.points.len() == 0 || !pt.eq(&contour.points[contour.points.len()-1])
                        {
                            contour.points.push(pt);
                        }
                        
                    }
                    // close contour if need be
                    if contour.points[0] != contour.points[contour.points.len()-1]
                    {
                        contour.points.push(contour.points[0].clone());
                    }
                    results.push(contour);
                    break;
                }                 
            }                
        }
        return results;
    }
}

impl<T> Contour for T where T: Triangulation + ?Sized, <T as Triangulation>::Vertex: GetValue{}

#[derive(Copy, Clone)]
struct VertexData 
{
    value: f64,
    position: Point2<f64>,
}
impl HasPosition for VertexData
{
    type Scalar=f64;

    fn position(&self) -> Point2<Self::Scalar> {
        self.position
    }
}
impl GetValue for VertexData
{
    fn get_value(&self) -> f64 {
        self.value
    }
}
impl Into<f64> for VertexData
{
    fn into(self) -> f64 {
       return self.value;
    }
}
fn f(x:f64, y:f64) -> f64{
    if x*x+ (y-50.0)*(y-50.0) < 25.0*25.0
    {
        return 25.0-f64::sqrt(x*x+ (y-50.0)*(y-50.0));
    }
    else if (x - 125.0)*(x - 125.0)+ (y-50.0)*(y-50.0) < 10.0*10.0
    {
        return 10.0-f64::sqrt((x - 125.0)*(x - 125.0)+ (y-50.0)*(y-50.0));
    }
    else {
        return 0.0
    }
}

#[test]
fn contour_test01()
{
    //let mut t:DelaunayTriangulation<VertexData> = DelaunayTriangulation::new();

    let mut points = Vec::new();
    points.reserve_exact(200*200);
    for i in 0..200
    {
        for j in (0..400).step_by(2)
        {
            let v = VertexData{position:Point2::new(0.5*j as f64, 0.5*i as f64), value: f(0.5*j as f64, 0.5*i as f64)};
            points.push(v);
            
        }
    }
    let t:DelaunayTriangulation<VertexData> = DelaunayTriangulation::<VertexData>::bulk_load(points).unwrap().0;
    let levels = vec!(1.0,5.0,10.0,20.0);
    for level in levels
    {
        let contours = t.contour(level, false);    
        for line in &contours
        {
            print!("\nLINESTRING(");
            for i in 0..line.points.len()
            {
                let delim = if i == line.points.len() - 1 {")"} else {", "};
                let point = line.points[i];
                print!("{} {}{delim}", point.x, point.y);
            }    
            print!("|{level}");
        }           
    }
}
    
