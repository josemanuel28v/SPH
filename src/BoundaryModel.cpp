#include "BoundaryModel.h" 

void BoundaryModel::init(std::vector<Vector3r> & points)
{
    resizeBoundary(points.size());

    for (unsigned int i = 0; i < points.size(); ++i)
        r[i] = points[i];
}

void BoundaryModel::resizeBoundary(const unsigned int size)
{
    r.resize(size);
}

void BoundaryModel::cleanBoundary()
{
    r.clear();
}
