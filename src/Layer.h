// Copyright 2022 Eric Fichter
#ifndef LAYER_H
#define LAYER_H

#include "headers.h"

//! Describes a geometric layer, e.g. located behind a space boundary.
//! The layer is linked to an ifc product.
//! Currently it is not used to represent a material layer but the ifc product itself.
//! It is defined by two distances, e.g. to the boundary
//! and a rank, that allows to prioritize over colliding layers.
//! It is assumed that layer faces are parallel to the boundary.
class Layer {

public:
    //! Constructor setting all member variables.
    Layer(Product *_product, double _dmin, double _dmax, double _rank);

    //! Returns begin of layer.
    double Dmin() const;

    //! Returns end of layer.
    double Dmax() const;

    //! Returns product.
    Product *RelProduct();

    //! Returns rank of layer.
    double Rank() const;

    //! Sets end of layer.
    void SetDmax(double _dmax);

private:
    //! Minimum distance, e.g. from space boundary to this.
    double dmin;

    //! Maximum distance, e.g. from space boundary to this.
    double dmax;

    //! Related product starting at dmin and ending at dmax.
    Product *product;

    //! Rank is currently calculated by negative of layer thickness. The thinner the layer, the higher the rank.
    //! So in case of collision of products, the smaller layer is kept.
    double rank;
};

#endif //LAYER_H