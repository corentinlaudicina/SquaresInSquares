#pragma once

#include "geometry/polygon.hpp"

#include <algorithm>  // std::min / std::max
#include <vector>

namespace pack {

// Implementation of a cell list.
//
// The square box [0, L]^2 is cut into nc x nc equal cells. Each particle is
// assigned the cell containing its centre. If every cell is at least as
// wide as the interaction range ("reach"), then two particles that interact
// must be in the same cell or in touching cells, so only a 3x3 block of
// cells needs to be searched.
//
// Particles in a cell are kept in a doubly linked list, stored in flat
// arrays (indices instead of pointers):
//
//   head[c]    -> first particle in cell c, or -1 if the cell is empty
//   next[i]    -> particle after i in its cell's list, or -1 at the end
//   prev[i]    -> particle before i, or -1 if i is the first (the head)
//   cell_of[i] -> the cell particle i is currently filed under
//
// Example: cell 5 holds particles 3, 7, 2 (in that order):
//   head[5]=3;  next[3]=7, next[7]=2, next[2]=-1;  prev[3]=-1, prev[7]=3, prev[2]=7
struct CellList {
    double side = 1.0;  // width of one cell (must be >= reach after build)
    int nc = 1;         // cells per side; total cells = nc * nc

    std::vector<int> head;     // per CELL: first particle, -1 if empty
    std::vector<int> next, prev, cell_of;  // per PARTICLE (see above)

    // Sets up the grid and files every particle into its cell. L = box side length, reach = max distance at which two particles interact.
    void build(const std::vector<Particle>& particles, double L, double reach);

    // Returns the cell index a centre at (x, y) belongs to.
    // Cells are numbered row by row: index = row * nc + column.
    int bin(double x, double y) const;

    // Call after particle i's centre changed. Re-files it only if it
    // crossed into a different cell.
    void move(int i, double x, double y);

    // Calls f(j) for every particle j in `cell` and its up-to-8 neighbouring
    // cells. F is any callable taking an int (lambda, function, functor). It is a template so the compiler can inline the call.
    template <class F>
    void for_each_near(int cell, F&& f) const {
        // Convert the 1D cell index back into (column, row).
        const int cx = cell % nc, cy = cell / nc;

        for (int y = std::max(cy - 1, 0); y <= std::min(cy + 1, nc - 1); ++y)
            for (int x = std::max(cx - 1, 0); x <= std::min(cx + 1, nc - 1); ++x)
                // Walk that cell's linked list: head -> next -> next ... -> -1
                for (int j = head[y * nc + x]; j >= 0; j = next[j])
                    f(j);
    }

    // Calls f(i, j) once for every pair of particles that are close enough
    // to possibly interact. 
    template <class F>
    void for_each_close_pair(F&& f) const {
        const int n = static_cast<int>(cell_of.size());
        for (int i = 0; i < n; ++i) {
            for_each_near(cell_of[i], [&](int j) {
                if (j > i) f(i, j);  // keep each unordered pair only once
            });
        }
    }

private:
    void link(int i, int cell);   // insert particle i at the front of a cell's list
    void unlink(int i);           // remove particle i from its current cell's list
};

}  // namespace pack