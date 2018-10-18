//  Copyright (C) 2002 Regents of the University of Michigan, portions used with
// permission
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//===================================================
//$Id$
//===================================================

#ifndef MDARRAY
#define MDARRAY

#include "stdio.h"

template <class T> class MDArray {
protected:
  T* data;
  int nDim; 
  long int nTotal; 
  long int *nSize_D;
  long int *nConvert_D;
  static const int nDimMax = 5; 

public:
  MDArray() {
    data = NULL;
    nDim = 0; 
    nSize_D = NULL; 
    nConvert_D = NULL; 
  };

  //===================================================
  ~MDArray() {
    if (data != NULL){
      delete[] data;
      delete[] nSize_D; 
      delete[] nConvert_D; 
    }
  };

  //===================================================
  MDArray(int nDimIn, long int n0, long int n1 = -1, long int n2 = -1, long int n3 = -1, long int n4 = -1){
    //init(int nDimIn, long int n0, long int n1, long int n2, long int n3, long int n4);
  };

  //===================================================
  void clear(){
    if (data != NULL){
      delete[] data;
      delete[] nSize_D; 
      delete[] nConvert_D; 
    }
    data = NULL;
    nDim = 0; 
    nSize_D = NULL; 
    nConvert_D = NULL;     
  }
  void init(int nDimIn, long int n0, long int n1 = -1, long int n2 = -1, long int n3 = -1, long int n4 = -1){
    nDim = nDimIn;

    // Check the dimension and size here!!!!!!!!!!!!!


    nSize_D = new long int[nDimMax];
    nConvert_D = new long int[nDimMax];
    for(int iDim = 0; iDim < nDimMax; iDim++){
      nSize_D[iDim] = 0; 
      nConvert_D[iDim] = 0; 
    }
        
    nSize_D[0] = n0; 
    if(n1 != -1) nSize_D[1] = n1;
    if(n2 != -1) nSize_D[2] = n2;
    if(n3 != -1) nSize_D[3] = n3;
    if(n4 != -1) nSize_D[4] = n4;


    nTotal = nSize_D[nDim-1]; 
    nConvert_D[nDim-1] = 1; 
    for(int iDim = nDim-2; iDim>=0; iDim--){
      nTotal *= nSize_D[iDim];
      nConvert_D[iDim] = nConvert_D[iDim+1]*nSize_D[iDim+1];
    }
      
    data = new T[nTotal];
    for (int i = 0; i<nTotal; i++) data[i] = 0; 
  };

  //===================================================
  T operator()(long int i0, long int i1=0, long int i2=0, long int i3=0,
               long int i4=0) const {
    return data[i0*nConvert_D[0] + 
		i1*nConvert_D[1] + 
		i2*nConvert_D[2] + 
		i3*nConvert_D[3] + 
		i4*nConvert_D[4]];
  };

  //===================================================
  T & operator()(long int i0, long int i1=0, long int i2=0, long int i3=0,
               long int i4=0){

    return data[i0*nConvert_D[0] + 
		i1*nConvert_D[1] + 
		i2*nConvert_D[2] + 
		i3*nConvert_D[3] + 
		i4*nConvert_D[4]];
  };

  void print(){
    for(int i=0; i<nTotal; i++){
      printf("idx = %i, value = %f\n", i, data[i]);
    }

  }

  //===================================================
};

#endif
