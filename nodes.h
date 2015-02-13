#ifndef NODE_H
#define NODE_H

#include <gmTrianglesystemModule>

class Nodes
{
public:
    Nodes();
    Nodes(GMlib::TSVertex<float> &vt);

    void setZ(float z);

    bool isThis(GMlib::TSVertex<float> *vt);

    GMlib::TSEdge<float>* getNeighbour(Nodes& nodes);
    GMlib::Array<GMlib::TSTriangle<float>* > getTriangles();


private:
    GMlib::TSVertex<float>* _vt;
};

#endif // NODE_H
