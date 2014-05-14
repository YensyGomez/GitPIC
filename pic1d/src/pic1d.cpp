#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <time.h>
#include <fstream>
#include <vector>
#include <fstream>
#include <fftw.h>


using namespace std;

//funciones
double distribution(double vb);
void Density(vector<double>r, vector<double>& n);
void fft_forward (vector<double>f, vector<double>&Fr,vector<double>&Fi);
void fft_backward (vector<double> Fr, vector<double> Fi,vector<double>& f);
void Electric (vector<double> phi, vector<double>& E);
void Poisson1D (vector<double>& u, vector<double> v, double kappa);
void Load (vector<double> r, vector<double> v, vector<double>& y);
void UnLoad (vector<double> y, vector<double>& r, vector<double>& v);
void rk4_fixed (double& x, vector<double>& y, void (*rhs_eval)(double, vector<double>, vector<double>&), double h);
void rhs_eval (double t, vector<double> y, vector<double>& dydt);
void Output (char* fn1, char* fn2, double t, vector<double> r, vector<double> v);



//parametros
double L; int N, C;


int main()
{
  // Parametros
  L =100.0;            // dominio de la solucion 0 <= x <= L (en longitudes de debye)
  N =100000;            // Numero de particulas
  C = 1000;            // Numero de celdas
  double vb = 3.0;    // velocidad rayo promedio
  double dt=0.1;    // delta tiempo (en frecuencias inversas del plasma)
  double tmax=1000;  // cantidad de iteraciones
  int skip = int (tmax / dt) / 10; //saltos del algoritmo para reportar datos

  ofstream vel;
  vel.open("datav.txt");

  ofstream pos;
  pos.open("datap.txt");

  ofstream graph;
    pos.open("graph.txt");



  vector<double> r, v, n(C); //r: posicion de las particulas, v: velocidad de particulas n: densidad de particulas por celda

  double t = 0.;
  int seed = time (NULL); srand (seed);
  for (int i = 0; i < N; i++)
  {
	  r.push_back(L*double (rand ()) / double (RAND_MAX));    //inicializando la posicion aleatoria
      v.push_back(distribution(vb));                          //inicializa la velocidad con una distribucion maxwelliana
  }



  for(int i=0; i<N;i++) //se imprimen los datos iniciales en un archivo
  {
	  vel<<v[i]<<endl;
	  pos<<r[i]<<endl;
	  graph<<r[i]<<" "<<v[i]<<endl;

  }
  vel.close();
  pos.close();

  char* phase[11]; char* data[11]; //archivos para almacenar los datos de salida
    phase[0] = "phase0.txt";phase[1] = "phase1.out";phase[2] = "phase2.out";
    phase[3] = "phase3.out";phase[4] = "phase4.out";phase[5] = "phase5.out";
    phase[6] = "phase6.out";phase[7] = "phase7.out";phase[8] = "phase8.out";
    phase[9] = "phase9.out";phase[10] = "phase10.out";data[0] = "data0.out";
    data[1] = "data1.out"; data[2] = "data2.out"; data[3] = "data3.out";
    data[4] = "data4.out"; data[5] = "data5.out"; data[6] = "data6.out";
    data[7] = "data7.out"; data[8] = "data8.out"; data[9] = "data9.out";
    data[10] = "data10.out";


    Output (phase[0], data[0], t, r, v); //inicializacion del algoritmo
    /*
     * La funcion Output calcula la densidad de las part�culas para cada celda
     * luego se calcula la ecuaci�n de poisson, la cual realiza una FFT que recibe el potencial
     * electrostatico. las ecuaciones estan explicadas en el paquete de fotocopias, en la parte donde dice
     * "soluci�n de la ecuacion de poisson" (normalizada)

     * */

  // iterando la solucion
      vector<double> y(2*N);
      Load (r, v, y);    //carga los valores de la posicion y la velocidad en el vector "y"
      for (int k = 1; k <= 10; k++)
        {
          for (int kk = 0; kk < skip; kk++)
            {
               // Take time-step
        	   rk4_fixed(t, y, rhs_eval, dt);

               // Make sure all coordinates in range 0 to L.
               for (int i = 0; i < N; i++)
                 {
                   if (y[i] < 0.) y[i] += L;
                   if (y[i] > L) y[i] -= L;
                 }

               //printf ("t = %11.4e\n", t);
            }
          //printf ("Plot %3d\n", k);

          // Output data
          UnLoad (y, r, v);
          Output(phase[k], data[k], t, r, v);
        }




      return 0;




}

double distribution (double vb)     //generador de distribuci�n maxwelliana para la velocidad
{
  // inicializa el generador aleatorio
  static int flag = 0;
  if (flag == 0)
    {
      int seed = time (NULL);
      srand (seed);
      flag = 1;
    }

  // Genera un valor random v
  double fmax = 0.5 * (1. + exp (-2. * vb * vb));
  double vmin = - 5. * vb;
  double vmax = + 5. * vb;
  double v = vmin + (vmax - vmin) * double (rand ()) / double (RAND_MAX);

  // Acceptar y reinyectar particulas
  double f = 0.5 * (exp (-(v - vb) * (v - vb) / 2.) +
		    exp (-(v + vb) * (v + vb) / 2.));
  double x = fmax * double (rand ()) / double (RAND_MAX);
  if (x > f) return distribution (vb);
  else
  {
	  //cout<<f<<endl;
	  return v;

  }

}
//c�lculo de la densidad de particulas en cada celda de la malla.
//si la part�cula esta en el borde del espacio de fase, coloca en la celda cero la densidad restante.

void Density (vector<double> r, vector<double>& n)
{
  // Initialize
  double dx = L / double (C);
  for(int x=0; x<C; x++)
  {
	  n[x]=0.0;
  }

  // Evaluar el numero de densidad
  for (int i = 0; i < N; i++)
    {
      int j = int (r[i] / dx);  //para saber en cual celda queda la particula (toma la parte entera)
      double y = r[i] / dx - double (j); // la posicion exacta de la particula dentro de la celda
      n[j] += (1. - y) / dx; //se le carga el valor a la celda de la diferencia y el
      if (j+1 == C) n[0] += y / dx;
      else n[j+1] += y / dx;
    }

  for(int x=0; x<C; x++)
    {
  	  //cout<<n[x]<<endl;
    }


}

// se usa la libreria fftw para calcular la transformada rapida de Fourier
// los vectores de entrada y salida son de tama�o C
// se calcula la fft del vector f en los vectores Fr y Fi

void fft_forward (vector<double>f, vector<double>&Fr, vector<double>&Fi)
{
  fftw_complex ff[C], FF[C]; //vectores complejos, donde ff[i].im es la parte imaginaria y ff[i].re es la parte real

  // se inician los valores de la parte real
  for (int j = 0; j < C; j++)
    {
      c_re (ff[j]) = f[j];
      c_im (ff[j]) = 0.;
    }

  //calculando el fft en 1 dimension
  fftw_plan p = fftw_create_plan (C, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_one (p, ff, FF);
  fftw_destroy_plan (p); //liberando memoria

  //copiando los resultados a la salida
  for (int j = 0; j < C; j++)
    {
      Fr[j] = c_re(FF[j]);
      Fi[j] = c_im(FF[j]);
    }

  // Normalizando la salida
  for (int j = 0; j < C; j++)
      {
        Fr[j]/=double (C);
        Fi[j]/=double (C);
      }

}

void fft_backward (vector<double> Fr, vector<double> Fi, vector<double>& f)
{
  fftw_complex ff[C], FF[C];

  //inicializacion
  for (int j = 0; j < C; j++)
    {
	  c_re(FF[j]) = Fr[j];
	  c_im(FF[j]) = Fi[j];
    }

  //calculando la fft inversa
  fftw_plan p = fftw_create_plan (C, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_one (p, FF, ff);
  fftw_destroy_plan (p); //liberando memoria

  //copiando los resultados a la salida
  for (int j = 0; j < C; j++)
      f[j] = c_re(ff[j]);
}

//The following routine solves Poisson's equation in 1-D to find the instantaneous electric potential on a uniform grid.

// Solves 1-d Poisson equation:
//    d^u / dx^2 = v   for  0 <= x <= L
// Periodic boundary conditions:
//    u(x + L) = u(x),  v(x + L) = v(x)
// Arrays u and v assumed to be of length J.
// Now, jth grid point corresponds to
//    x_j = j dx  for j = 0,J-1
// where dx = L / J.
// Also,
//    kappa = 2 pi / L

void Poisson1D (vector<double>& u, vector<double> v, double kappa)
{

  vector<double> Vr(C), Vi(C), Ur(C), Ui(C);

  for (int i=0; i<Vr.size(); i++)
    {
  	  Vr[i]=0.;
  	  Vi[i]=0.;
  	  Ur[i]=0.;
  	  Ui[i]=0.;
    }

  // Fourier transform source term
  fft_forward (v, Vr, Vi);


  // calcula fft para u
  Ur[0] = Ui[0] = 0.;
  for (int j=1;j<=C/2;j++)
    {
      Ur[j] = - Vr[j] / double (j * j) / kappa / kappa;
      Ui[j] = - Vi[j] / double (j * j) / kappa / kappa;
    }
  for (int j = C/2;j<C;j++)
    {
      Ur[j] = Ur[C-j];
      Ui[j] = - Ui[C-j];
    }

  // fft inversa para hallar u
  fft_backward (Ur, Ui, u);


}

// Calculate electric field from potential

void Electric (vector<double> phi, vector<double>& E)
{
  double dx = L / double (C);


  for (int j = 1; j < C-1; j++)
  {
	  E[j] = (phi[j-1] - phi[j+1]) / 2. / dx;
  }
  E[0] = (phi[C-1] - phi[1]) / (2. * dx);
  E[C-1] = (phi[C-2] - phi[0]) /  (2. * dx);
}
// Ecuaciones de movimiento.
// Electron equations of motion:
//    y(0:N-1)  = r_i
//    y(N:2N-1) = dr_i/dt

void rhs_eval (double t, vector<double> y, vector<double>& dydt)
{
  // Declare local arrays
  vector<double> r(N), v(N), rdot(N), vdot(N), r0(N);
  vector<double> ne(C), rho(C), phi(C), E(C);

  // Unload data from y
  //UnLoad (y, r, v);

  // Make sure all coordinates in range 0 to L
  r0 = r;
  for (int i = 0; i < N; i++)
    {
      if (r0[i] < 0.) r0[i] += L;
      if (r0[i] > L) r0[i] -= L;
    }

  // Calculate electron number density
  Density (r0, ne);



  // Solve Poisson's equation
  double n0 = double (N) / L;
  for (int j = 0; j < C; j++)
    rho[j] = ne[j] / n0 - 1.;
  double kappa = 2. * M_PI / L;
  Poisson1D (phi, rho, kappa);

  //for(int i=0;i<phi.size();i++) cout<<phi[i]<<endl;
  // Calculate electric field
  Electric (phi, E);

  // Equations of motion
  for (int i = 0; i < N; i++)
    {
      double dx = L / double (C);
      int j = int (r0[i] / dx);
      double y = r0[i] / dx - double (j);

      double Efield;
      if (j+1 == C)
         Efield = E[j] * (1. - y) + E[0] * y;
      else
         Efield = E[j] * (1. - y) + E[j+1] * y;

      rdot[i] = v[i];
      vdot[i] = - Efield;
    }

  // Load data into dydt
  Load (rdot, vdot, dydt);
}

void Load (vector<double> r, vector<double> v, vector<double>& y)
{
  for (int i = 0; i < N; i++)
    {
      y[i] = r[i];
      y[N+i] = v[i];
    }
}

// Unload particle coordinates from solution vector

void UnLoad (vector<double> y, vector<double>& r, vector<double>& v)
{
  for (int i = 0; i < N; i++)
    {
      r[i] = y[i];
      v[i] = y[N+i];
    }
}

void rk4_fixed (double& x, vector<double>& y,
                void (*rhs_eval)(double, vector<double>, vector<double>&),
                double h)
{
  // Array y assumed to be of extent n, where n is no. of coupled
  // equations
  int n = y.size();

  // Declare local arrays
  vector<double> k1(n), k2(n), k3(n), k4(n), f(n), dydx(n);

  // Zeroth intermediate step
  (*rhs_eval) (x, y, dydx);
  for (int j = 0; j < n; j++)
    {
      k1[j] = h * dydx[j];
      f[j] = y[j] + k1[j] / 2.;
    }

  // First intermediate step
  (*rhs_eval) (x + h / 2., f, dydx);
  for (int j = 0; j < n; j++)
    {
      k2[j] = h * dydx[j];
      f[j] = y[j] + k2[j] / 2.;
    }

  // Second intermediate step
  (*rhs_eval) (x + h / 2., f, dydx);
  for (int j = 0; j < n; j++)
    {
      k3[j] = h * dydx[j];
      f[j] = y[j] + k3[j];
    }

  // Third intermediate step
  (*rhs_eval) (x + h, f, dydx);
  for (int j = 0; j < n; j++)
    {
      k4[j] = h * dydx[j];
    }

  // Actual step
  for (int j = 0; j < n; j++)
    {
      y[j] += k1[j] / 6. + k2[j] / 3. + k3[j] / 3. + k4[j] / 6.;
    }
  x += h;

  return;
}

void Output (char* fn1, char* fn2, double t,
	     vector<double> r, vector<double> v)
{
  // Write phase-space data
  FILE* file = fopen (fn1, "w");
  for (int i = 0; i < N; i++)
    fprintf (file, "%e %e\n", r[i], v[i]);
  fclose (file);

  // Write electric field data
  vector<double> ne(C), n(C), phi(C), E(C);
  for (int i=0; i<n.size(); i++)
  {
	  ne[i]=0.;
	  n[i]=0.;
	  phi[i]=0.;
	  E[i]=0.;
  }

  Density (r, ne);


  for (int j = 0; j < C; j++)
  {
	  n[j] = double (C) * ne[j] / double (N) - 1.;
  }
  double kappa = 2. * M_PI / L;

  Poisson1D (phi, n, kappa);
  Electric (phi, E);

  file = fopen (fn2, "w");
  for (int j = 0; j < C; j++)
    {
      double x = double (j) * L / double (C);
      fprintf (file, "%e %e %e %e\n", x, ne[j], n[j], E[j]);
    }
  double x = L;
  fprintf (file, "%e %e %e %e\n", x, ne[0], n[0], E[0]);
  fclose (file);
}


