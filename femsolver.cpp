#include "femsolver.h"

#include "drum.h"

/**
 * @brief FEMSolver::FEMSolver
 * @param n Triangles on innermost circle
 * @param m Number of circles
 * @param radius Radius of circle
 * @param scene
 */
FEMSolver::FEMSolver(int n, int m, float radius, std::shared_ptr<GMlib::Scene> scene)
{
    GMlib::ArrayLX<GMlib::TSVertex<float> > temp;
    makeRegular(temp, n, m, radius);
    _scene = scene;

    _theDrum = new Drum(temp);
    _theDrum->setMaxForce(0.1);
    _theDrum->prepareComputation();
    _theDrum->computeNormals();
    _theDrum->replot();
    _scene->insert(_theDrum);
}

void FEMSolver::makeRegular(GMlib::ArrayLX<GMlib::TSVertex<float> > &array, int n, int m, float radius)
{
    //array.setSize(1 + n * (m * (m + 1) ) / 2);
    //array[0].setPos(GMlib::Point<float,2>(0,0));
    array += GMlib::TSVertex<float>(GMlib::Point<float,2>(0,0) );

    for(int i = 0; i < m; i++)
    {
        int inner = n * (i + 1);

        GMlib::Angle a(M_2PI / (inner - 1) );
        GMlib::SqMatrix<float,2> matrix(a, GMlib::Vector<float,2>(1,0), GMlib::Vector<float,2>(0,1) );
        GMlib::Point<float,2> point(radius/m, 0);

        array += GMlib::TSVertex<float>(point);
        for(int j = 1; j < inner; j++)
        {
            point = matrix * point;
            array += GMlib::TSVertex<float>(point);
        }
    }
}

void FEMSolver::makeRandom(GMlib::ArrayLX<GMlib::TSVertex<float> > &array, int n, float radius)
{

}
