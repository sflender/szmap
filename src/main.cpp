#include "healpix_map.h"
#include "healpix_map_fitsio.h"
#include "fitsio.h"
#include "config.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cosmo.h"
#include "main.h"

using namespace std;

int main(int argc, char *argv[]){

    if(argc<2){
        cout<<"---usage:\n";
        cout<<"---./szmap config_file"; 
        cout<<"------\n";
        return -1;
    }

    // --- paramters from input file
    ConfigFile cfg(argv[1]);
	int order = cfg.getValueOfKey<int>("order");
	int order_subpix = cfg.getValueOfKey<int>("order_subpix");
	bool convert_map_to_ring = cfg.getValueOfKey<bool>("convert_map_to_ring");
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
    float max_int_radius = cfg.getValueOfKey<float>("max_int_radius"); // where to truncate the profile (in units of R200c)
    float trunc_radius = cfg.getValueOfKey<float>("trunc_radius");
    float Mmin = cfg.getValueOfKey<float>("Mmin");
    float Mmax = cfg.getValueOfKey<float>("Mmax");
    float zmin = cfg.getValueOfKey<float>("zmin");
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
	//int rank = MPI::COMM_WORLD.Get_rank();
	//int nranks = MPI::COMM_WORLD.Get_size();
	int rank = MPI::COMM_WORLD.Get_rank();
	int nranks = MPI::COMM_WORLD.Get_size();


	if (rank==0) cout<<nranks<<" MPI ranks "<<endl;

	// --- Healpix parameters for the SZ map
	int pix_mult = pow(4,order_subpix-order); // pixel multiplier (how many subpixels are in a main pixel)

	int nside = pow(2,order);
	long int nside_subpix = pow(2,order_subpix);
	int npix = 12*pow(nside,2);
	long int npix_subpix = 12*pow(nside_subpix,2);

	Healpix_Base HPB = Healpix_Base(order,NEST);
	Healpix_Base2 HPB_subpix = Healpix_Base2(order_subpix,NEST);

	if(rank==0){
		cout<<"npix = "<<npix<<endl;
		cout<<"sub-pixelized npix = "<<npix_subpix<<endl;
		cout<<"pixel multiplier = "<<pix_mult<<endl;
	}

	// --- pixel arrays to store the map
	float *local_pixel_array = new float[npix]; fill(local_pixel_array,local_pixel_array+npix,0.0);
	float *global_pixel_array    = new float[npix]; fill(global_pixel_array,global_pixel_array+npix,0.0);

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

	cosmo cosm_model(H0, Omega_M, Omega_b, Omega_k, w_DE); // construct the cosm. model (needed for ang diam distance)

   	vector<int> sel;
   	int nsel=0;
   	if(rank==0){cout<<"selecting halos\n";}
   	for (int i=0; i<N_clus; i++){
   		if( (M200red[i]>Mmin and M200red[i]<Mmax and redshift[i]>zmin and redshift[i]<zmax) ){
   			sel.push_back(i);
   			nsel++;
   		}
   	}

   	if(rank==0){cout<<nsel<<" selected halos. now painting in SZ profiles\n";}



   	// --- now loop over all the halos
   	int completed = 0;
	for (int j=rank; j<nsel; j+=nranks){ // distribute the halos evenly over all MPI ranks
		int i = sel[j];

		float d_obs_halo = cosm_model.ang_diam(redshift[i]); // physical distance in Mpc (=ang.diameter distance) from observer to halo
		float rhocrit = cosm_model.calc_rho_crit(redshift[i]);
        float hubblez = cosm_model.hubblez(redshift[i]);
        float Delta_vir = cosm_model.Delta_vir(redshift[i]);
        float Omega_Mz = cosm_model.Omega_Mz(redshift[i]);

		// --- convert from ra/dec to Healpix coord's:
		float theta = (90.0 - dec[i])*deg2rad;
		float phi = ra[i] * deg2rad;

        printf("%f %f", theta,phi);
		pointing ptg(theta,phi);
		ptg.normalize();
		vec3 vec_to_cluster = ptg.to_vec3(); // vector pointing to the cluster

        float M200 = M200red[i]/h; //M200 in Msol
        float R200 = pow(M200/((4.0/3.0)*M_PI*200*rhocrit), 1.0/3.0);
        float theta200 = R200/d_obs_halo;
        float theta_max = trunc_radius * theta200;

   		vector<int> *listpix = new vector<int>; //list of pixel indices in disc generated by query_disc
	    HPB.query_disc_inclusive(ptg, theta_max, *listpix);
		int npix_in_disc = listpix->size();

		for (int p=0; p<npix_in_disc; p++){ // loop over all pixels in the disc
			
			for (long int q = (long int)pix_mult* (long int)(*listpix)[p]; q < (long int)pix_mult*( (long int)(*listpix)[p]+1); q++){ // loop over all subpixels inside the main pixel
				
                printf("q %d\n", q);
				pointing ptg_pix = HPB_subpix.pix2ang(q);
				ptg_pix.normalize();
				float alpha = acos( dotprod(vec_to_cluster,ptg_pix.to_vec3() ) ); // angular separation in rad		
                //printf("alpha %f \n", alpha);
				float d_halo_pixel = alpha * d_obs_halo; // physical distance from halo center to pixel p (Mpc)
				//float ysz = return_ysz_Battaglia(max_int_radius, d_halo_pixel/R200, M200, R200, redshift[i], rhocrit, Omega_b/Omega_M);
                //float ysz = return_ysz_Beta_virial(M200, redshift[i], h, Omega_M, Omega_Mz, Delta_vir, d_obs_halo);
                float ysz = return_ysz_Beta(max_int_radius, d_halo_pixel, M200, redshift[i], h, Omega_M, Omega_Mz, Delta_vir, d_obs_halo, rhocrit);

				float tsz = ysz * tsz_norm;
				local_pixel_array[ (*listpix)[p]] += tsz / pix_mult; //increment tsz at pixel p by the tsz fraction generated by subpixel q
			}
		}

		delete listpix;

		if (completed%1==0 and rank==0) cout<<"rank "<<rank<<" completed tSZ profile "<<completed<<" of "<<nsel/nranks<<"\n";
		completed++;
	}

	if(rank==0){cout<<"reducing tsz map\n";}
	MPI_Reduce(local_pixel_array, global_pixel_array, npix, MPI_FLOAT, MPI_SUM, 0, comm);

	delete [] ra, dec, M200red, redshift, local_pixel_array;

	if (rank==0 and write_output_file==1){
		cout << "...writing tsz output fits file "<<output_tsz_name<<endl;

		Healpix_Map<float> *HPmap = new Healpix_Map<float>;
		if (convert_map_to_ring==0){
			HPmap->SetNside(nside,NEST);
			for (int p=0;p<npix;p++){
				(*HPmap)[p] = global_pixel_array[p];
			}
		}

/*
		else if (convert_map_to_ring==1){
			HPmap->SetNside(nside,RING);
			for (int p=0;p<npix;p++){
				(*HPmap)[p] = global_pixel_array[HPB.nest2ring(p)];
			}
		}
*/
		
		else if (convert_map_to_ring==1){
			HPmap->SetNside(nside,RING);
			for (int pnest=0;pnest<npix;pnest++){
				int pring = HPB.nest2ring(pnest);
				(*HPmap)[pring] = global_pixel_array[pnest];
			}
		}


		write_Healpix_map_to_fits(fitsout_tsz, *HPmap, PLANCK_FLOAT32);
		
		delete HPmap;

	}
	
	if(rank==0) cout<<"...all done\n";

 	delete [] global_pixel_array;

	MPI_Finalize();

	return 0;
}
