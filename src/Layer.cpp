// Copyright 2022 Eric Fichter
#include "Layer.h"

Layer::Layer(Product *_product, double _dmin, double _dmax, double _rank) : product(_product), dmin(_dmin), dmax(_dmax), rank(_rank) {}

double Layer::Dmin() const { return dmin; }

double Layer::Dmax() const { return dmax; }

Product *Layer::RelProduct() { return product; }

double Layer::Rank() const { return rank; }

void Layer::SetDmax(double _dmax) { dmax = _dmax; }