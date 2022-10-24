#ifndef slic3r_NonplanarFacet_hpp_
#define slic3r_NonplanarFacet_hpp_

#include "libslic3r.h"
#include "Geometry.hpp"

namespace Slic3r {

typedef struct {
  float x;
  float y;
  float z;
} facet_vertex;

typedef struct {
  facet_vertex    max;
  facet_vertex    min;
} facet_stats;

class NonplanarFacet
{
    public:
    facet_vertex vertex[3];
    facet_vertex normal;
    float U[3];
    float V[3];
    float N[3];
    float a[3];
    float b[3];
    float dot_product;
    float maga;
    float magb;
    float theta;
    int neighbor[3];
    facet_stats stats;
    bool marked = false;

    NonplanarFacet() {};
    ~NonplanarFacet() {};
    void calculate_stats();
    void calculate_theta(float t);
    void translate(float x, float y, float z);
    void scale(float versor[3]);
    float calculate_surface_area();

};
};

#endif
