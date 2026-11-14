#pragma once

// Initial conditions: the overlap-free configurations a run starts a System from, before any
// mover has touched it. These are not movers -- a move changes a configuration that already
// exists; grid_init() and diagonal_init() build one from nothing but a shape and a count.

#include "system.hpp"

namespace pack {

System grid_init(const Shape& shape, int n);

System diagonal_init(const Shape& shape, int n);

}  // namespace pack
