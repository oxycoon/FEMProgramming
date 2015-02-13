#include "nodes.h"

Nodes::Nodes()
{
}

Nodes::Nodes(GMlib::TSVertex<float> &vt)
{
    _vt = &vt;
}

void Nodes::setZ(float z)
{
    _vt->setZ(z);
}

bool Nodes::isThis(GMlib::TSVertex<float> *vt)
{
    return _vt == vt;
}

GMlib::TSEdge<float> *Nodes::getNeighbour(Nodes &nodes)
{
    GMlib::ArrayT<GMlib::TSEdge<float>*> edge = _vt->getEdges();

    for(int i = 0; i < edge.size(); i++)
    {
        if(nodes.isThis(edge[i]->getOtherVertex(*_vt) ) )
        {
            return edge[i];
        }
    }
    return NULL;
}

GMlib::Array<GMlib::TSTriangle<float> *> Nodes::getTriangles()
{
    return _vt->getTriangles();
}
