#include "healpix_map.h"
#include "healpix_map_fitsio.h"
#include "fitsio.h"
#include "config.h"
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cosmo.h"
#include "main.h"

using namespace std;

int main(){

    if(argc<2){
        cout<<"---usage:\n";
        cout<<"---./makemap config_file"; 
        cout<<"------\n";
        return -1;
    }

    // --- paramters from input file
    ConfigFile cfg(argv[1]);
	int order = cfg.getValueOfKey<int>("order");
	int order_subpix = cfg.getValueOfKey<int>("order_subpix");
	bool write_output_file = cfg.getValueOfKey<bool>("write_output_file");
	string output_tsz_name = cfg.getValueOfKey<string>("output_tsz_name");
	string cluster_catalog  = cfg.getValueOfKey<string>("cluster_catalog");
	float nu = cfg.getValueOfKey<float>("nu");

    //---cosmo inputs
    float h = cfg.getValueOfKey<float>("h");
    float Omega_M = cfg.getValueOfKey<float>("Omega_M");
    float Omega_b = cfg.getValueOfKey<float>("Omega_b");
    float Omega_k = cfg.getValueOfKey<float>("Omega_k");
    float w_DE = cfg.getValueOfKey<float>("w_DE");

    //---other
    float max_int_radius = cfg.getValueOfKey<float>("max_int_radius");
    float Mmin = cfg.getValueOfKey<float>("Mmin");
    float Mmax = cfg.getValueOfKey<float>("Mmax");
    float zmin = xfg.getValueOfKey<float>("zmin");
    float zmax = cfg.getValueOfKey<float>("zmax");

	// --- hard-coded constants + derived params
	float T_CMB_mikroK = 2.725*1e6;
	float rad2deg = 180.0/M_PI;
	float deg2rad = M_PI/180.0;
	float rad2arcmin = 180.0*60.0/M_PI;
	float arcmin2rad = 1.0/rad2arcmin;
	float conc_norm=1.0; //additonal parameter for cM relation
	float conc_mass_norm=1.0; //additonal parameter for cM relation

	float H0=100.0*h;
	float xx=nu/56.68;
	float g_nu=xx/tanh(xx/2.0)-4.0;
	float tsz_norm = g_nu * T_CMB_mikroK;

	// --- set MPI initialization
	MPI_Init(&argc,&argv);
	MPI_Comm comm = MPI_COMM_WORLD;
	int rank = MPI::COMM_WORLD.Get_rank();
	int nranks = MPI::COMM_WORLD.Get_size();

	// --- Healpix parameters for the tSZ map
	int npix = 12*nside*nside;
	Healpix_Map<float> *HPmap = new Healpix_Map<float>;
	HPmap->SetNside(nside,RING);
	HPmap->fill(0.0);

	// --- additional pixel arrays (those are needed because cannot call MPI_REDUCE on Healpix maps)
	float *local_pixel_array_tsz = new float[npix]; fill(local_pixel_array_tsz,local_pixel_array_tsz+npix,0);
	float *global_pixel_array    = new float[npix]; fill(global_pixel_array,global_pixel_array+npix,0);

	// --- prepare the output
	fitshandle fitsout_ksz, fitsout_tsz;
	if (rank==0 and write_output_file==1){
		fitsout_tsz.create(output_tsz_name);
	}

	// --- read in the cluster catalog
	// some definitions for cfitsio
	fitsfile *fptr;
	int status=0;
	int nhdus, current_hdu, hdutype;
	int ncols=0;
	long int N_clus;
	long int repeat, width;
	int casesen = 0; // case sensitivity
	int colnum;
	int typecode;
	fits_open_file(&fptr, cluster_catalog.c_str(), READONLY, &status);
   	if (status)
   	{
      printf("unable to read files\n");
      return -1;
   	}

	fits_movabs_hdu(&fptr[0], 2, &hdutype, &status); //move to HDU 2 (that's where the data sits)
	fits_get_num_rows(&fptr[0], &N_clus, &status); // get the total number of clusters in the catalog

	// read in the data
	float *ra = new float[N_clus];
	fits_get_colnum(&fptr[0], casesen, (char*)"RA", &colnum, &status);
	fits_get_coltype(&fptr[0], colnum, &typecode, &repeat, &width, &status);
	fits_read_col(&fptr[0], typecode, colnum, 1, 1, N_clus, NULL, ra, NULL, &status);

	float *dec = new float[N_clus];
	fits_get_colnum(&fptr[0], casesen, (char*)"DEC", &colnum, &status);
	fits_get_coltype(&fptr[0], colnum, &typecode, &repeat, &width, &status);
	fits_read_col(&fptr[0], typecode, colnum, 1, 1, N_clus, NULL, dec, NULL, &status);

	float *M200red = new float[N_clus];
	fits_get_colnum(&fptr[0], casesen, (char*)"M200RED", &colnum, &status);
	fits_get_coltype(&fptr[0], colnum, &typecode, &repeat, &width, &status);
	fits_read_col(&fptr[0], typecode, colnum, 1, 1, N_clus, NULL, M200red, NULL, &status);

	float *redshift = new float[N_clus];
	fits_get_colnum(&fptr[0], casesen, (char*)"REDSHIFT", &colnum, &status);
	fits_get_coltype(&fptr[0], colnum, &typecode, &repeat, &width, &status);
	fits_read_col(&fptr[0], typecode, colnum, 1, 1, N_clus, NULL, redshift, NULL, &status);

	return 0;
}