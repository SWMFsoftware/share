#ifndef _GRIDINFO_H_
#define _GRIDINFO_H_

#include <cmath>
#include <iostream>

extern "C" {
void reset_grid_info(int* gridInfo, int nx, int ny, int nz);

void set_point_status(int* gridInfo, int nx, int ny, int nz, int i, int j,
                      int k, int val);

void get_point_status(int* gridInfo, int nx, int ny, int nz, int i, int j,
                      int k, int* Val);
}

int get_point_status(int* gridInfo, int nx, int ny, int nz, int i, int j,
                     int k);

class GridInfo {
public:
  constexpr static int iPicOn_ = 1, iPicOff_ = 0;

private:
  int nx, ny, nz, nPatch, nInt;
  int* data;
  bool isNewGrid;

public:
  GridInfo() {
    data = nullptr;
    isNewGrid = false;
  }

  ~GridInfo() {
    if (data != nullptr) {
      delete[] data;
    }
  }

  bool is_grid_new() const { return isNewGrid; }

  int get_patch_size(int iDim) {
    int size = 1;
    if (iDim == 0)
      size = nPatch;
    if (iDim == 1 && ny > 1)
      size = nPatch;
    if (iDim == 2 && nz > 1)
      size = nPatch;
    return size;
  }

  void init(int nxIn, int nyIn, int nzIn, int nPatchIn) {
    nx = nxIn;
    ny = nyIn;
    nz = nzIn;
    nPatch = nPatchIn;

    nx = nx / nPatch;
    if (ny > 1)
      ny = ny / nPatch;
    if (nz > 1)
      nz = nz / nPatch;

    int nInt = ceil(nx * ny * nz / 8.0 / sizeof(int));
    data = new int[nInt];

    for (int i = 0; i < nInt; i++) {
      data[i] = 0;
    }
    isNewGrid = true;
  }

  // void reset() { reset_grid_info(data, nx, ny, nz); }

  // void set_status(int i, int j, int k, int val) {
  //   set_point_status(data, nx, ny, nz, i, j, k, val);
  // }

  void set_status(int* statusIn) {
    if (statusIn == nullptr) {
      turn_on_all_cells();
      return;
    }
    const int nInt = ceil(nx * ny * nz / 8.0 / sizeof(int));
    isNewGrid = false;
    for (int i = 0; i < nInt; i++) {
      const int tmp = statusIn[i];
      if (data[i] != tmp)
        isNewGrid = true;
      data[i] = tmp;
    }
  }

  void turn_on_all_cells() {
    int iOn = 1;
    for (int i = 0; i < nx; i++)
      for (int j = 0; j < ny; j++)
        for (int k = 0; k < nz; k++) {
          set_point_status(data, nx, ny, nz, i, j, k, iOn);
        }
    return;
  }

  int get_status(int i, int j, int k) {
    // i, j, k: FLEKS cell index. NOT patch index.
    int iPatch = i / nPatch, jPatch = j, kPatch = k;
    if (ny > 1)
      jPatch /= nPatch;
    if (nz > 1)
      kPatch /= nPatch;
    return get_point_status(data, nx, ny, nz, iPatch, jPatch, kPatch);
  }
};

#endif
