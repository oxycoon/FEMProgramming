#ifndef DRUM_H
#define DRUM_H

#include <gmTrianglesystemModule>

#include "nodes.h"

class Drum : public GMlib::TriangleFacets<float>
{
public:
    Drum(int n, int m, float radius);
    Drum(const GMlib::ArrayLX<GMlib::TSVertex<float> > &array);

    void setForce(float force);
    void setMaxForce(float maxForce);

    float getForce() const;
    float getMaxForce() const;

    void prepareComputation();

protected:
    void localSimulate(double dt);

private:
    float _force;
    float _maxForce;

    int _direction; //Force direction

    GMlib::ArrayLX<Nodes> _nodes;
    GMlib::DMatrix<float> _a; //Stiffness matrix
    GMlib::DVector<float> _b; //load vector

    GMlib::Vector<GMlib::Vector<float,2>, 3> findVectors(GMlib::TSEdge<float>* edge);
    GMlib::Vector<GMlib::Vector<float,2>, 3> findVectors(GMlib::TSTriangle<float>* triangle, Nodes* nodes);

    void updateHeights(float force);

};

#endif // DRUM_H
