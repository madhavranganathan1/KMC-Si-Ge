#include "variables.h"

int update(AtomX,AtomY)
{
  int AtomZ=Height[AtomX][AtomY];
  double toss = gsl_rng_uniform_pos(gnalea_r4);
            
  int Bond=0;
                 
        int TYPExyz=displacementType[AtomX][AtomY][AtomZ].Type;
                
                             displacementType[AtomX][AtomY][AtomZ].Type = 0 ;
                             Height[AtomX][AtomY] -= 1;

                if (toss < 0.25) {
                             displacementType[(AtomX -1 + Nx)% Nx][AtomY][(Height[(AtomX -1 + Nx)% Nx][AtomY]) + 1].Type = TYPExyz;
                             Height[(AtomX -1 + Nx)% Nx][AtomY] += 1 ;
                             Bond=1+Totalbond((AtomX-1+Nx)%Nx,AtomY);
                                 }

                else if (toss < 0.5) {
                             displacementType[(AtomX+1) % Nx][AtomY][(Height[(AtomX+1) % Nx][AtomY]) + 1].Type = TYPExyz;
                             Height[(AtomX +1)% Nx][AtomY] += 1;
                             Bond=1+Totalbond((AtomX+1)%Nx,AtomY);
                                      }

                else if (toss < 0.75) {
                            displacementType[AtomX][(AtomY-1+Ny) % Ny][(Height[AtomX][(AtomY-1+Ny) % Ny]) + 1].Type = TYPExyz;
                            Height[AtomX][(AtomY - 1 + Ny)%Ny] += 1 ;
                            Bond=1+Totalbond(AtomX,(AtomY-1+Ny)%Ny);
                                       }

                else if (toss < 1.0001) {
                             displacementType[AtomX][(AtomY + 1) % Ny][(Height[AtomX][(AtomY + 1) % Ny]) + 1].Type = TYPExyz;
                             Height[AtomX][(AtomY +1)% Ny] += 1;
                             Bond = 1+Totalbond(AtomX,(AtomY+1)%Ny);
                                        }


return Bond;
}
