#ifndef FEMSOLVER_H
#define FEMSOLVER_H

#include <gmTrianglesystemModule>

#include "nodes.h"

#include <memory>




class FEMSolver : public GMlib::TriangleFacets<float>
{
public:
    FEMSolver(int d = 0);
    FEMSolver(const GMlib::ArrayLX<GMlib::TSVertex<float> > &array);
    //FEMSolver(int n, int m, float radius, std::shared_ptr<GMlib::Scene> scene);

    void makeRegular(int n, int m, float radius);
    void makeRandom(int n, float radius);

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

#endif // FEMSOLVER_H
