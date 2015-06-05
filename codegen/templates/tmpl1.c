// initialisation
for(int i=0; i<=Nx; i++)
    u[i][0] = I(x[i])
// separate calculation for first time step
for(int i=1; i<Nx; i++)
    u[i][1] = u[i][0]
        -0.5*pow(C,2)*(u[i+1][0] - 2*u[i][0] + u[i-1][0])
// main loop
for(int n=1; n<Nt; n++){
    // boundary conditions
    u[n+1][0] = 0; u[n+1][Nx] = 0;
    for(int i=1; i<Nx; i++)
    // update mesh points at time = n+1
        ${loopbody}
}
