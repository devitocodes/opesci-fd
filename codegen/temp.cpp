    #pragma omp single
    {
    printf("%d\n", _ti);
    for(int x=0;x<dimx;x++){
      for(int y=0;y<dimy;y++){
        printf("%.5f ", Txx[t][x][y]);
      }
      printf("\n");
    }
    }