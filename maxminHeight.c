int maxminHeight()
{
  int In,J,H1D[Nx*Ny];

  for(In=0;In<Nx;In++){
    for(J=0;J<Ny;J++){
       H1D[In*Nx+J]=Height[In][J];}}
 
  int max = H1D[0];
  int min = H1D[0];

  for (In=0; In<Nx*Ny; In++)
        {
          if (H1D[In] > max)
                {
                  max = H1D[In];
                }
          else if (H1D[In] < min)
                {
                  min = H1D[In];
                }
        }
  if(rank==1){printf ("Max H %d and Min H: %d \n", max,min);}

return 0;
}
