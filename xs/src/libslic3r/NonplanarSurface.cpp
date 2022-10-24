#include "NonplanarSurface.hpp"
#include <unordered_set>
#include <iostream>
#include <fstream>


namespace Slic3r {

NonplanarSurface::NonplanarSurface(std::map<int, NonplanarFacet> &_mesh)
{
    this->mesh = _mesh;
    this->calculate_stats();
}

bool
NonplanarSurface::operator==(const NonplanarSurface& other) const {
   return (stats.min.x == other.stats.min.x &&
           stats.min.y == other.stats.min.y &&
           stats.min.z == other.stats.min.z &&
           stats.max.x == other.stats.max.x &&
           stats.max.y == other.stats.max.y &&
           stats.max.z == other.stats.max.z);
}

void
NonplanarSurface::calculate_stats()
{
    //calculate min and max values
    this->stats.min.x = 10000000;
    this->stats.min.y = 10000000;
    this->stats.min.z = 10000000;
    this->stats.max.x = -10000000;
    this->stats.max.y = -10000000;
    this->stats.max.z = -10000000;
    for(auto& facet : this->mesh) {
        this->stats.min.x = std::min(this->stats.min.x, facet.second.stats.min.x);
        this->stats.min.y = std::min(this->stats.min.y, facet.second.stats.min.y);
        this->stats.min.z = std::min(this->stats.min.z, facet.second.stats.min.z);
        this->stats.max.x = std::max(this->stats.max.x, facet.second.stats.max.x);
        this->stats.max.y = std::max(this->stats.max.y, facet.second.stats.max.y);
        this->stats.max.z = std::max(this->stats.max.z, facet.second.stats.max.z);
    }
}

void
NonplanarSurface::translate(float x, float y, float z)
{
    //translate all facets
    for(auto& facet : this->mesh) {
        facet.second.translate(x, y, z);
    }

    //translate min and max values
    this->stats.min.x += x;
    this->stats.min.y += y;
    this->stats.min.z += z;
    this->stats.max.x += x;
    this->stats.max.y += y;
    this->stats.max.z += z;
}

void
NonplanarSurface::scale(float factor)
{
    float versor[3];
    versor[0] = factor;
    versor[1] = factor;
    versor[2] = factor;
    this->scale(versor);
}

void
NonplanarSurface::scale(float versor[3])
{
    //scale all facets
    for(auto& facet : this->mesh) {
        facet.second.scale(versor);
    }

    //scale min and max values
    this->stats.min.x *= versor[0];
    this->stats.min.y *= versor[1];
    this->stats.min.z *= versor[2];
    this->stats.max.x *= versor[0];
    this->stats.max.y *= versor[1];
    this->stats.max.z *= versor[2];
}

void
NonplanarSurface::rotate_z(float angle) {
  double radian_angle = (angle / 180.0) * 3.1415927;
  double c = cos(radian_angle);
  double s = sin(radian_angle);

  for(auto& facet : this->mesh) {
    for(int j = 0; j < 3; j++) {
        double xold = facet.second.vertex[j].x;
        double yold = facet.second.vertex[j].y;
        facet.second.vertex[j].x = c * xold - s * yold;
        facet.second.vertex[j].y = s * xold + c * yold;
    }
    facet.second.calculate_stats();
  }

  this->calculate_stats();
}


   
void
NonplanarSurface::debug_output()
{
    std::cout << "Facets(" << this->mesh.size() << "): (min:X:" << this->stats.min.x << " Y:" << this->stats.min.y << " Z:" << this->stats.min.z <<
                           " max:X:" << this->stats.max.x << " Y:" << this->stats.max.y << " Z:" << this->stats.max.z << ")" <<
                           "Height " << this->stats.max.z - this->stats.min.z << std::endl;
    for(auto& facet : this->mesh) {
        // Every facet has a normal which represents the orientation of the facet, z-comp. can be used to calc. facet angle
        std::cout << "triangle: (" << facet.first << ")(" << facet.second.marked << ") ";
        std::cout << " (" << (180*std::acos(facet.second.normal.z))/3.14159265 << "°)";

        std::cout << " | V0:";
        std::cout << " X:"<< facet.second.vertex[0].x;
        std::cout << " Y:"<< facet.second.vertex[0].y;
        std::cout << " Z:"<< facet.second.vertex[0].z;

        std::cout << " | V1:";
        std::cout << " X:"<< facet.second.vertex[1].x;
        std::cout << " Y:"<< facet.second.vertex[1].y;
        std::cout << " Z:"<< facet.second.vertex[1].z;

        std::cout << " | V2:";
        std::cout << " X:"<< facet.second.vertex[2].x;
        std::cout << " Y:"<< facet.second.vertex[2].y;
        std::cout << " Z:"<< facet.second.vertex[2].z;

        std::cout << " | Normal:";
        std::cout << " X:"<< facet.second.normal.x;
        std::cout << " Y:"<< facet.second.normal.y;
        std::cout << " Z:"<< facet.second.normal.z;



        //TODO check if neighbors exist
        // stl_neighbors* neighbors = mesh.stl.neighbors_start + facet.first;
        std::cout << " | Neighbors:";
        std::cout << " 0:"<< facet.second.neighbor[0];
        std::cout << " 1:"<< facet.second.neighbor[1];
        std::cout << " 2:"<< facet.second.neighbor[2];
        std::cout << std::endl;
    }
}
// calculate the normal of each facet (in case the normal function does not work, this is definetly the correct formula)
void 
NonplanarSurface::calculate_normal(){

    for(auto& facet : this->mesh){
        facet.second.U[0] = (facet.second.vertex[1].x - facet.second.vertex[0].x);
        facet.second.U[1] = (facet.second.vertex[1].y - facet.second.vertex[0].y);
        facet.second.U[2] = (facet.second.vertex[1].z - facet.second.vertex[0].z);
        facet.second.V[0] = (facet.second.vertex[2].x - facet.second.vertex[0].x);
        facet.second.V[1] = (facet.second.vertex[2].y - facet.second.vertex[0].y);
        facet.second.V[2] = (facet.second.vertex[2].z - facet.second.vertex[0].z);

        ///facet.second.N[0] = (facet.second.vertex[0].y*facet.second.vertex[1].z) - (facet.second.vertex[0].z*facet.second.vertex[1].y);
        ///facet.second.N[1] = (facet.second.vertex[0].z*facet.second.vertex[1].x) - (facet.second.vertex[0].x*facet.second.vertex[1].z);
        ///facet.second.N[2] = (facet.second.vertex[0].x*facet.second.vertex[1].y) - (facet.second.vertex[0].y*facet.second.vertex[1].x);
        facet.second.N[0] = (facet.second.U[1]*facet.second.V[2]) - (facet.second.U[2]*facet.second.V[1]);
        facet.second.N[1] = (facet.second.U[2]*facet.second.V[0]) - (facet.second.U[0]*facet.second.V[2]);
        facet.second.N[2] = (facet.second.U[0]*facet.second.V[1]) - (facet.second.U[1]*facet.second.V[0]);
        std::ofstream myfile;
        std::stringstream stream;
        
        myfile.open ("/home/shantiverschoor/workspaces/Slic3r/normal.txt", std::ios_base::app);
        stream <<std::fixed <<std::setprecision(3)<<"x "<<facet.second.N[0]<<"y "<< facet.second.N[1]<<"z " <<facet.second.N[2]<<"\n";
        stream <<std::fixed <<std::setprecision(3)<<"vertex 0: "<<facet.second.vertex[0].x<<" vertex 1: "<<facet.second.vertex[1].x<<"vertex 2: "<<facet.second.vertex[2].x<<"\n";
        stream <<std::fixed <<std::setprecision(3)<<"vertex 0: "<<facet.second.vertex[0].y<<" vertex 1: "<<facet.second.vertex[1].y<<"vertex 2: "<<facet.second.vertex[2].y<<"\n";
        stream <<std::fixed <<std::setprecision(3)<<"vertex 0: "<<facet.second.vertex[0].z<<" vertex 1: "<<facet.second.vertex[1].z<<"vertex 2: "<<facet.second.vertex[2].z<<"\n";
        myfile <<stream.str();
        myfile.close();
    }

}
void
NonplanarSurface::dot_product(){
    for(auto& facet : this->mesh){
        facet.second.a[0] = facet.second.N[0];
        facet.second.a[1] = facet.second.N[1];
        facet.second.a[2] = facet.second.N[2];
        facet.second.b[0] = 0;
        facet.second.b[1] = 0;
        facet.second.b[2] = 1;
        facet.second.dot_product = facet.second.a[0]*facet.second.b[0] +facet.second.a[1]*facet.second.b[1]+facet.second.a[2]*facet.second.b[2];
        std::ofstream myfile;
        std::stringstream stream;
        myfile.open ("/home/shantiverschoor/workspaces/Slic3r/dot.txt", std::ios_base::app);
        stream <<std::fixed <<std::setprecision(3)<< facet.second.dot_product<<"\n";
        myfile <<stream.str();
        myfile.close();
  
    
}
    }
}

void
NonplanarSurface::mag(){
    for(auto& facet : this->mesh){
        facet.second.a[0] = facet.second.N[0];
        facet.second.a[1] = facet.second.N[1];
        facet.second.a[2] = facet.second.N[2];
        facet.second.b[0] = 0;
        facet.second.b[1] = 0;
        facet.second.b[2] = 1;
        facet.second.maga = std::sqrt(facet.second.a[0]*facet.second.a[0]+facet.second.a[1]*facet.second.a[1]+facet.second.a[2]*facet.second.a[2]);
        facet.second.magb = std::sqrt(facet.second.b[0]*facet.second.b[0]+facet.second.b[1]*facet.second.b[1]+facet.second.b[2]*facet.second.b[2]);
        std::ofstream myfile;
        std::stringstream stream;
        myfile.open ("/home/shantiverschoor/workspaces/Slic3r/mag.txt", std::ios_base::app);
        stream <<std::fixed <<std::setprecision(3)<<"a="<<facet.second.maga <<" b="<<facet.second.magb<<"\n";
        myfile <<stream.str();
        myfile.close();
    
}
}

//here lies my plans to make a theta vector :/
///Theta
//NonplanarSurface::theta(){
   // for(auto& facet : this->mesh){
   //     theta.push_back(facet.second.theta);
  //  }
//}

void
NonplanarSurface::calculate_theta(){
    // the first facet gets a theta value (for now, can be looped over all facets)
    //if surface is flat (test cheeseslice) use the first facet's theta as surface theta
    //for more complicated surfaces, you need to do something with multiple thetas (tbc)
 auto facet = this->mesh.begin();
    /// update: to make this return theta it needs to be set in the function calculate_theta and called in PrintObject. Otherwise the function does not know how to get correct normal, dot etc.
    facet->second.U[0] = (facet->second.vertex[1].x - facet->second.vertex[0].x);
    facet->second.U[1] = (facet->second.vertex[1].y - facet->second.vertex[0].y);
    facet->second.U[2] = (facet->second.vertex[1].z - facet->second.vertex[0].z);
    facet->second.V[0] = (facet->second.vertex[2].x - facet->second.vertex[0].x);
    facet->second.V[1] = (facet->second.vertex[2].y - facet->second.vertex[0].y);
    facet->second.V[2] = (facet->second.vertex[2].z - facet->second.vertex[0].z);
    facet->second.N[0] = (facet->second.U[1]*facet->second.V[2]) - (facet->second.U[2]*facet->second.V[1]);
    facet->second.N[1] = (facet->second.U[2]*facet->second.V[0]) - (facet->second.U[0]*facet->second.V[2]);
    facet->second.N[2] = (facet->second.U[0]*facet->second.V[1]) - (facet->second.U[1]*facet->second.V[0]);
    facet->second.a[0] = facet->second.N[0];
    facet->second.a[1] = facet->second.N[1];
    facet->second.a[2] = facet->second.N[2];
    facet->second.b[0] = 0;
    facet->second.b[1] = 0;
    facet->second.b[2] = 1;
    
     
    facet->second.dot_product = facet->second.a[0]*facet->second.b[0] +facet->second.a[1]*facet->second.b[1]+facet->second.a[2]*facet->second.b[2];
    facet->second.maga = std::sqrt(facet->second.a[0]*facet->second.a[0]+facet->second.a[1]*facet->second.a[1]+facet->second.a[2]*facet->second.a[2]);
    facet->second.magb = std::sqrt(facet->second.b[0]*facet->second.b[0]+facet->second.b[1]*facet->second.b[1]+facet->second.b[2]*facet->second.b[2]);

    this->theta = 180* std::acos((facet->second.dot_product/(facet->second.maga*facet->second.magb)))/3.14159265;
    
    std::ofstream myfile;
    std::stringstream stream;
    myfile.open ("/home/shantiverschoor/workspaces/Slic3r/thetas.txt", std::ios_base::app);
    stream <<std::fixed <<std::setprecision(3)<< "dot/mag= "<<facet->second.dot_product/(facet->second.maga*facet->second.magb)<<" a*b = "<<facet->second.maga*facet->second.magb<<"theta=  "<<this->theta<<"\n";
    myfile <<stream.str();
    myfile.close();
    
}


NonplanarSurfaces
NonplanarSurface::group_surfaces()
{
    std::pair<int, NonplanarFacet> begin = *this->mesh.begin();
    this->mark_neighbor_surfaces(begin.first);
    NonplanarSurface newSurface;
    for (std::map<int, NonplanarFacet>::iterator it=this->mesh.begin(); it!=this->mesh.end();) {
        if((*it).second.marked == false) {
            newSurface.mesh[(*it).first] = (*it).second;
            it = this->mesh.erase(it);
        }
        else {
            ++it;
        }
    }
    this->calculate_stats();
    if (newSurface.mesh.size() == 0) {
        //return only this
        NonplanarSurfaces nonplanar_surfaces;
        nonplanar_surfaces.push_back(*this);
        return nonplanar_surfaces;
    }
    else {
        //return union of this and recursion
        NonplanarSurfaces nonplanar_surfaces = newSurface.group_surfaces();
        nonplanar_surfaces.push_back(*this);
        return nonplanar_surfaces;
    }
}

void
NonplanarSurface::mark_neighbor_surfaces(int id)
{
    //if already marked return
    if(this->mesh.find(id) == this->mesh.end() || this->mesh[id].marked == true) return;

    std::unordered_set<int> todo;
    todo.insert(id);
    while (!todo.empty()) {
        int id = *todo.begin();
        todo.erase(todo.begin());
        //mark this facet
        this->mesh[id].marked = true;
        //mark all neighbors
        for(int j = 0; j < 3; j++) {
            int neighbor_id = this->mesh[id].neighbor[j];
            if(this->mesh.find(neighbor_id) != this->mesh.end() && !this->mesh[neighbor_id].marked) {
                todo.insert(neighbor_id);
            }
        }
    }
}

bool
NonplanarSurface::check_max_printing_height(float height)
{
    if ((this->stats.max.z - this->stats.min.z) > height ) {
        std::cout << "Surface removed: printheight too heigh (" << (this->stats.max.z - this->stats.min.z) << " mm)" << '\n';
        return true;
    } else {
        return false;
    }
}

bool
NonplanarSurface::check_surface_area(float minimal_area)
{
    //calculate surface area of nonplanar surface.
    float area = 0.0f;
    for(auto& facet : this->mesh) {
        area += facet.second.calculate_surface_area();
    }
    if (area < minimal_area) {
        std::cout << "Surface removed: area too small (" << area << " mm²)" << '\n';
        return true;
    } else {
        return false;
    }
}

void
NonplanarSurface::check_printable_surfaces(float max_angle)
{
    //TODO do something
}

/* this will return scaled ExPolygons */
ExPolygons
NonplanarSurface::horizontal_projection() const
{
    Polygons pp;
    pp.reserve(this->mesh.size());
    for(auto& facet : this->mesh) {
        Polygon p;
        p.points.resize(3);
        p.points[0] = Point(facet.second.vertex[0].x / SCALING_FACTOR, facet.second.vertex[0].y / SCALING_FACTOR);
        p.points[1] = Point(facet.second.vertex[1].x / SCALING_FACTOR, facet.second.vertex[1].y / SCALING_FACTOR);
        p.points[2] = Point(facet.second.vertex[2].x / SCALING_FACTOR, facet.second.vertex[2].y / SCALING_FACTOR);
        p.make_counter_clockwise();  // do this after scaling, as winding order might change while doing that
        pp.push_back(p);
    }

    // the offset factor was tuned using groovemount.stl
    return union_ex(offset(pp, 0.01 / SCALING_FACTOR), true);
}

