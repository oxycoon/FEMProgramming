#ifndef FEMSOLVER_H
#define FEMSOLVER_H

#include <gmTrianglesystemModule>

#include "nodes.h"

#include <memory>

class Drum;


class FEMSolver
{
public:
    FEMSolver(int n, int m, float radius, std::shared_ptr<GMlib::Scene> scene);

    void makeRegular(GMlib::ArrayLX<GMlib::TSVertex<float> > &array, int n, int m, float radius);
    void makeRandom(GMlib::ArrayLX<GMlib::TSVertex<float> > &array, int n, float radius);

private:
    Drum* _theDrum;

    std::shared_ptr<GMlib::Scene> _scene;
};

#endif // FEMSOLVER_H
