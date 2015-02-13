#include "femsolver.h"

#include "drum.h"

/**
 * @brief FEMSolver::FEMSolver
 * @param n Triangles on innermost circle
 * @param m Number of circles
 * @param radius Radius of circle
 * @param scene
 */
/*FEMSolver::FEMSolver(int n, int m, float radius, std::shared_ptr<GMlib::Scene> scene)
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
}*/

FEMSolver::FEMSolver(int d) : TriangleFacets<float>(d)
{

}

FEMSolver::FEMSolver(const GMlib::ArrayLX<GMlib::TSVertex<float> > &array): TriangleFacets<float>(array)
{

}

void FEMSolver::makeRegular(int n, int m, float radius)
{
    GMlib::ArrayLX<GMlib::TSVertex<float> > array;
    array.setSize(1 + n * (m * (m + 1) ) / 2);
    array[0].setPos(GMlib::Point<float,2>(0,0));

    int elements = 1;
    for(int i = 0; i < m; i++)
    {
        int inner = n * (i + 1);
        for(int j = 0; j < inner; j++)
        {
            GMlib::Angle a(j*M_2PI / (inner) );
            GMlib::SqMatrix<float,2> matrix(a, GMlib::Vector<float,2>(1,0), GMlib::Vector<float,2>(0,1) );
            array[elements++].setPos(matrix * GMlib::Point<float,2>(radius * (i+1)/m, 0));
        }
    }
    this->insertAlways(array);
}

void FEMSolver::makeRandom(int n, float radius)
{

}

void FEMSolver::setForce(float force)
{
     _force = force;
}

void FEMSolver::setMaxForce(float maxForce)
{
    _maxForce = maxForce;
}

float FEMSolver::getForce() const
{
    return _force;
}

float FEMSolver::getMaxForce() const
{
    return _maxForce;
}

void FEMSolver::prepareComputation()
{
    _direction = 1;

    GMlib::Vector<GMlib::Vector<float,2>, 3> d;
    GMlib::Array<GMlib::TSTriangle<float>*> triangles;

    float dd; //diagonal for d vector

    //triangle 1
    float area1;
    float dh1;
    float h1;

    //triangle 2
    float area2;
    float dh2;
    float h2;


    //Perform delaunay
    this->triangulateDelaunay();
    for(int i = 0; i < this->size(); i++)
    {
        if(!(*this)[i].boundary())
        {
            _nodes += Nodes((*this)[i]);
        }
    }

    //initialize stiffness matrix
    _a.setDim(_nodes.size(), _nodes.size());
    for(int i = 0; i < _nodes.size(); i++)
    {
        for(int j = 0; j < _nodes.size(); j++)
        {
            _a[i][j] = 0;
        }
    }

    //Create stiffness matrix
    for(int i = 0; i < _nodes.size(); i++)
    {
        for(int j = 0; j < _nodes.size(); j++)
        {
            GMlib::TSEdge<float> *currentEdge = _nodes[i].getNeighbour(_nodes[j]);

            if(currentEdge != NULL)
            {
                d = findVectors(currentEdge);
                dd = 1/(d[0] * d[0]);

                //triangle 1
                dh1 = dd * (d[1] * d[0]);
                area1 = std::abs(d[0] ^ d[1]);
                h1 = dd * area1;

                //triangle 2
                dh2 = dd * (d[2] * d[0]);
                area2 = std::abs(d[0] ^ d[2]);
                h2 = dd * area1;

                //Stiffness elements outside diagonal
                _a[i][j] = _a[j][i] = ( (dh1 * (1 - dh1) ) / h1 - dd) * area1 / 2 +
                                      ( (dh2 * (1 - dh2) ) / h2 - dd) * area2 / 2;
            }
        }
        //Compute diagonal for stiffness matrix
        triangles = _nodes[i].getTriangles();
        for(int j = 0; j < triangles.size(); j++)
        {
            d = findVectors(triangles[j], &_nodes[i]);
            _a[i][i] += d[2] * d[2] / (2 * std::abs(d[0] ^ d[1]) );
        }
    }
    _a.invert();

    std::cout << _a << std::endl;

    //Create load vectors
    _b.setDim(_nodes.size());
    for(int i = 0; i < _nodes.size(); i++)
    {
        GMlib::Array<GMlib::TSTriangle<float>* > triangles2 = _nodes[i].getTriangles();
        _b[i] = triangles2[0]->getArea2D() / 3;
        for(int j = 0; j < triangles2.size(); j++)
        {
            _b[i] += triangles2[j]->getArea2D() / 3;
        }
    }
}

void FEMSolver::localSimulate(double dt)
{
    _force += (dt / 8 ) * _direction;

    if(_force > _maxForce)
        _direction = -1;
    if(_force < -_maxForce)
        _direction = 1;

    updateHeights(_force);

    this->replot();
}

void FEMSolver::updateHeights(float force)
{
    GMlib::DVector<float> b = force * _b;
    GMlib::DVector<float> x = _a * b;

    for(int i = 0; i < _nodes.size(); i++)
    {
        _nodes[i].setZ(x[i]);
    }
}

GMlib::Vector<GMlib::Vector<float, 2>, 3> FEMSolver::findVectors(GMlib::TSTriangle<float> *triangle, Nodes *node)
{
    GMlib::Vector<GMlib::Vector<float,2>, 3> d; //return vector

    GMlib::Array<GMlib::TSVertex<float>* > v = triangle->getVertices();
    GMlib::Point<float,2> p0, p1, p2; // points of triangle

    //Control if points are in the correct order compared to the node.
    if(node->isThis(v[1])) //Check if the node is the 2nd point of the triangle
    {
        std::swap(v[0], v[1]);
        std::swap(v[1], v[2]);
    }

    if(node->isThis(v[2])) //Check if the node is the 3rd point of the triangle
    {
        std::swap(v[0], v[2]);
        std::swap(v[1], v[2]);
    }

    p0 = v[0]->getParameter();
    p1 = v[1]->getParameter();
    p2 = v[2]->getParameter();

    d[0] = p2 - p0;
    d[1] = p1 - p0;
    d[2] = p2 - p1;

    return d;
}

GMlib::Vector<GMlib::Vector<float, 2>, 3> FEMSolver::findVectors(GMlib::TSEdge<float> *edge)
{
    GMlib::Vector<GMlib::Vector<float,2>, 3> d; //return vector

    GMlib::Array<GMlib::TSTriangle<float>* > triangle = edge->getTriangle();
    GMlib::Array<GMlib::TSVertex<float>* > v1 = triangle[0]->getVertices();
    GMlib::Array<GMlib::TSVertex<float>* > v2 = triangle[1]->getVertices();

    GMlib::Point<float,2> p0, p1, p2, p3; // points of triangles

    p0 = edge->getFirstVertex()->getParameter();
    p1 = edge->getLastVertex()->getParameter();

    for(int i = 0; i < 3; i++)
    {
        if(v1[i] != edge->getFirstVertex() && v1[i] != edge->getLastVertex())
        {
            p2 = v1[i]->getParameter();
        }
        if(v2[i] != edge->getFirstVertex() && v2[i] != edge->getLastVertex())
        {
            p3 = v2[i]->getParameter();
        }
    }

    //Vector definitions for 2D:
    d[0] = p1 - p0; //d vector
    d[1] = p2 - p0; //a1 vector
    d[2] = p3 - p0; //a2 vector

    return d;
}
