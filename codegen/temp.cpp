
std::vector<float> u(2*dimx*dimy*dimz), v(2*dimx*dimy*dimz), w(2*dimx*dimy*dimz),
    txx(2*dimx*dimy*dimz), tyy(2*dimx*dimy*dimz), tzz(2*dimx*dimy*dimz),
    tyz(2*dimx*dimy*dimz), txz(2*dimx*dimy*dimz), txy(2*dimx*dimy*dimz);

    float (*U)[dimx][dimy][dimz] = (float (*)[dimx][dimy][dimz]) u.data();
    float (*V)[dimx][dimy][dimz] = (float (*)[dimx][dimy][dimz]) v.data();
    float (*W)[dimx][dimy][dimz] = (float (*)[dimx][dimy][dimz]) w.data();

    float (*Txx)[dimx][dimy][dimz] = (float (*)[dimx][dimy][dimz]) txx.data();
    float (*Tyy)[dimx][dimy][dimz] = (float (*)[dimx][dimy][dimz]) tyy.data();
    float (*Tzz)[dimx][dimy][dimz] = (float (*)[dimx][dimy][dimz]) tzz.data();
    float (*Tyz)[dimx][dimy][dimz] = (float (*)[dimx][dimy][dimz]) tyz.data();
    float (*Txz)[dimx][dimy][dimz] = (float (*)[dimx][dimy][dimz]) txz.data();
    float (*Txy)[dimx][dimy][dimz] = (float (*)[dimx][dimy][dimz]) txy.data(); 




printf("%d\n", _ti);
for(int i=2;i<dimx - 2;++i){
    for(int j=2;j<dimy - 2;++j){
      for(int k=2;k<dimz - 2;++k){
        printf("%.4f ",V[0][i][j][k]);
      }
      printf("\n");
    }
    printf("\n");
}

U[t1][x][1][z]=((  - dy*V[t1][x][1][z] + dy*V[t1][x][2][z] + dy*V[t1][x + 1][1][z] - dy*V[t1][x + 1][2][z])
U[t1][x][dimy - 2][z]=((  - dy*V[t1][x][dimy - 4][z] + dy*V[t1][x][dimy - 3][z] + dy*V[t1][x + 1][dimy - 4][z] - dy*V[t1][x + 1][dimy - 3][z])
