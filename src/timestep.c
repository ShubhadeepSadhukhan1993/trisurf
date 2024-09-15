/* vim: set ts=4 sts=4 sw=4 noet : */
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
//#include "io.h"
#include "general.h"
#include "timestep.h"
#include "vertexmove.h"
#include "bondflip.h"
#include "frame.h"
#include "io.h"
#include "stats.h"
#include "sh.h"
#include "shcomplex.h"
#include "vesicle.h"
#include<gsl/gsl_complex.h>
#include<gsl/gsl_complex_math.h>
#include<string.h>
#include <sys/stat.h>
#include "cluster.h"
#include "vertex.h"

#include "energy.h"
/*Shubhadeep*/

double error_check(ts_vesicle *vesicle){
	ts_double x1, y1, z1,  x2, y2, z2, vmag;
	x1=vesicle->xnorm*vesicle->vmag;
	y1=vesicle->ynorm*vesicle->vmag;
	z1=vesicle->znorm*vesicle->vmag;

	//printf("vmag2 %1.5e\n", vesicle->vmag);

	/*
	int ii, jj;
	x2=0;
	y2=0;
	z2=0;
	for (jj=0; jj<vesicle->vlist->n; jj++){
		if (vesicle->vlist->vtx[jj]->c>1e-15){
			x2+=vesicle->vlist->vtx[jj]->Factx;
			y2+=vesicle->vlist->vtx[jj]->Facty;
			z2+=vesicle->vlist->vtx[jj]->Factz;
		}
	}
	*/
	ts_double dir[3];
	if (vesicle->tape->inhibition_switch && vesicle->tape->number_of_vertices_with_c0!=0 && vesicle->tape->F!=0){
		calculate_vesicle_dir(vesicle, dir);
	}
	calculate_velocity_mag(vesicle,  &vmag);	

	x2=vmag*dir[0];
	y2=vmag*dir[1];
	z2=vmag*dir[2];
	//printf("old-> %1.2f %1.2f %1.2f\n",x1, y1, z1 );
	//printf("new-> %1.2f %1.2f %1.2f\n",x2, y2, z2 );
	return (pow(x1-x2, 2)+pow(y1-y2, 2)+pow(z1-z2, 2))/(pow(x1,2)+pow(y1,2)+pow(z1, 2));

}



ts_double find_max(double x[3]){
	ts_double xmax=-1e300;
	int i;
	for (i=0; i<3; i++){
		if (x[i]>=xmax){
			xmax=x[i];
		}
	}
	return xmax;
}


ts_double vector_error(ts_double x[3], double y[3]){
	ts_double ans1=0;
	ts_double ans2=0;
	int i;
	for (i=0; i<3; i++){
		ans1+=pow(x[i]-y[i],2);
		ans2+=pow(x[i],2);	
	}
	return ans1/ans2;

}

void matrix_multiplication(double A[3][3], ts_double x[3]){
	ts_double x_new[3];
	x_new[0]=x[0];
	x_new[1]=x[1];
	x_new[2]=x[2];

	x[0]=0;
	x[1]=0;
	x[2]=0;
	int i,j;
	for (i=0; i<3; i++){
		for (j=0; j<3; j++){
			x[i]+=A[i][j]*x_new[j];
		}
	}
}

void find_vesicle_dir(ts_vesicle *vesicle, double A[3][3], double x[3]){
	int i, j;
	x[0]=drand48()+1;
	x[1]=drand48()+1;
	x[2]=drand48()+1;
	double x_old[3];
	ts_double err, xmax;

	/*for (i=0; i<3; i++){
		for (j=0; j<3; j++){
			printf("%1.e\t",A[i][j]);
		}
		printf("\n");
	}*/
	while (1){
		for (i=0; i<3; i++){x_old[i]=x[i];}
		matrix_multiplication(A, x);
		
		xmax=find_max(x);
		x[0]/=xmax;
		x[1]/=xmax;
		x[2]/=xmax;
		//ts_double norm =0;
		//for (i=0; i<3; i++){norm+=pow(x[i],2);}
		//norm=sqrt(norm);
	
		//x[0]/=norm;
		//x[1]/=norm;
		//x[2]/=norm;
		//printf("%1.ef \t %1.e \t %1.e\n",x_old[0], x_old[1], x_old[2]);
		err=vector_error(x, x_old);
		//printf("%1.ef\n",err);
		if (err<1e-3){break;}
	}
	ts_double norm =0;
	for (i=0; i<3; i++){norm+=pow(x[i],2);}
	norm=sqrt(norm);
	
	x[0]/=norm;
	x[1]/=norm;
	x[2]/=norm;
	
}


void calculate_velocity_mag(ts_vesicle *vesicle, double *vmag){
	ts_double x0=0;
	ts_double y0=0;
	ts_double z0=0;
	int i;
	for (i=0; i<vesicle->vlist->n; i++){
		if (vesicle->vlist->vtx[i]->c>1e-10){
			x0+=vesicle->xnorm*vesicle->vlist->vtx[i]->Factx;
			y0+=vesicle->ynorm*vesicle->vlist->vtx[i]->Facty;
			z0+=vesicle->znorm*vesicle->vlist->vtx[i]->Factz;
		}
	}
	*vmag=(x0+y0+z0);
	if (*vmag<0){
		vesicle->xnorm*=-1;
		vesicle->ynorm*=-1;
		vesicle->znorm*=-1;
		*vmag*=-1;
	}

}


void calculate_vesicle_dir(ts_vesicle *vesicle, double dir[3]){
	
	ts_cluster_list * cstlist=init_cluster_list();
	clusterize_vesicle(vesicle, cstlist);
	
	ts_double *zlist=calloc(cstlist->n, sizeof(ts_double));
	ts_double *ylist=calloc(cstlist->n, sizeof(ts_double));
	ts_double *xlist=calloc(cstlist->n, sizeof(ts_double));

	
	int i;
	for ( i=0; i<cstlist->n; i++){
		cluster_mean(cstlist->cluster[i], vesicle);
		xlist[i]=cstlist->cluster[i]->x;
		ylist[i]=cstlist->cluster[i]->y;
		zlist[i]=cstlist->cluster[i]->z;
	}


	ts_double axx, axy, azx, ayy, ayz, azz;
	axx=0; axy=0; ayz=0; azx=0; ayy=0; azz=0;

	ts_double A[3][3];

	for (i=0; i<cstlist->n; i++){
		
		axx+=cstlist->cluster[i]->nvtx*(xlist[i]-vesicle->cm[0])*(xlist[i]-vesicle->cm[0]);
		ayy+=cstlist->cluster[i]->nvtx*(ylist[i]-vesicle->cm[1])*(ylist[i]-vesicle->cm[1]);
		azz+=cstlist->cluster[i]->nvtx*(zlist[i]-vesicle->cm[2])*(zlist[i]-vesicle->cm[2]);
		axy+=cstlist->cluster[i]->nvtx*(xlist[i]-vesicle->cm[0])*(ylist[i]-vesicle->cm[1]);
		ayz+=cstlist->cluster[i]->nvtx*(ylist[i]-vesicle->cm[1])*(zlist[i]-vesicle->cm[2]);
		azx+=cstlist->cluster[i]->nvtx*(zlist[i]-vesicle->cm[2])*(xlist[i]-vesicle->cm[0]);
	}
	
	A[0][0]=axx;
	A[1][1]=ayy;
	A[2][2]=azz;

	A[0][1]=axy;
	A[1][0]=axy;

	A[1][2]=ayz;
	A[2][1]=ayz;
	
	A[2][0]=azx;
	A[0][2]=azx;
	
	find_vesicle_dir(vesicle, A, dir);

	cluster_list_free(cstlist);
	free_vtx_cluster(vesicle);

}


void calculate_vesicle_dir1(ts_vesicle *vesicle, double dir[3]){
	ts_double *xlist=calloc(vesicle->tape->number_of_vertices_with_c0, sizeof(ts_double));
	ts_double *ylist=calloc(vesicle->tape->number_of_vertices_with_c0, sizeof(ts_double));
	ts_double *zlist=calloc(vesicle->tape->number_of_vertices_with_c0, sizeof(ts_double));
	ts_double *wlist=calloc(vesicle->tape->number_of_vertices_with_c0, sizeof(ts_double));
	int i, j;
	j=0;
	for ( i=0; i<vesicle->vlist->n; i++){
		if (vesicle->vlist->vtx[i]->c>1e-10){
			xlist[j]=vesicle->vlist->vtx[i]->x;
			ylist[j]=vesicle->vlist->vtx[i]->y;
			zlist[j]=vesicle->vlist->vtx[i]->z;
			wlist[j]=sqrt(pow(vesicle->vlist->vtx[i]->Factx,2)+pow(vesicle->vlist->vtx[i]->Facty,2)+pow(vesicle->vlist->vtx[i]->Factz,2));
			j+=1;
		}
	}

	ts_double axx, axy, azx, ayy, ayz, azz;
	axx=0; axy=0; ayz=0; azx=0; ayy=0; azz=0;

	ts_double A[3][3];

	for (i=0; i<vesicle->tape->number_of_vertices_with_c0; i++){
		axx+=wlist[i]*(xlist[i]-vesicle->cm[0])*(xlist[i]-vesicle->cm[0]);
		ayy+=wlist[i]*(ylist[i]-vesicle->cm[1])*(ylist[i]-vesicle->cm[1]);
		azz+=wlist[i]*(zlist[i]-vesicle->cm[2])*(zlist[i]-vesicle->cm[2]);
		axy+=wlist[i]*(xlist[i]-vesicle->cm[0])*(ylist[i]-vesicle->cm[1]);
		ayz+=wlist[i]*(ylist[i]-vesicle->cm[1])*(zlist[i]-vesicle->cm[2]);
		azx+=wlist[i]*(zlist[i]-vesicle->cm[2])*(xlist[i]-vesicle->cm[0]);
	}
	
	A[0][0]=axx;
	A[1][1]=ayy;
	A[2][2]=azz;

	A[0][1]=axy;
	A[1][0]=axy;

	A[1][2]=ayz;
	A[2][1]=ayz;
	
	A[2][0]=azx;
	A[0][2]=azx;
	find_vesicle_dir(vesicle, A, dir);

}


void calculate_vesicle_dir2(ts_vesicle *vesicle, double dir[3]){
	
	ts_cluster_list * cstlist=init_cluster_list();
	clusterize_vesicle(vesicle, cstlist);
	
	ts_double *zlist=calloc(cstlist->n, sizeof(ts_double));
	ts_double *ylist=calloc(cstlist->n, sizeof(ts_double));
	ts_double *xlist=calloc(cstlist->n, sizeof(ts_double));
	ts_double *wlist=calloc(cstlist->n, sizeof(ts_double));

	
	int i;
	for ( i=0; i<cstlist->n; i++){
		cluster_mean(cstlist->cluster[i], vesicle);
		xlist[i]=cstlist->cluster[i]->x;
		ylist[i]=cstlist->cluster[i]->y;
		zlist[i]=cstlist->cluster[i]->z;
		wlist[i]=cstlist->cluster[i]->nvtx;//sqrt(pow(cstlist->cluster[i]->Factx,2)+pow(cstlist->cluster[i]->Facty,2)+pow(cstlist->cluster[i]->Factz,2));
		//cstlist->cluster[i]->Factx*(cstlist->cluster[i]->x-vesicle->cm[0])+cstlist->cluster[i]->Facty*(cstlist->cluster[i]->y-vesicle->cm[1]);//+cstlist->cluster[i]->Factz*(cstlist->cluster[i]->z-vesicle->cm[2]);
		
		//printf("%e\n", wlist[i]);
		//sqrt(pow(cstlist->cluster[i]->Factx,2)+pow(cstlist->cluster[i]->Facty,2)+pow(cstlist->cluster[i]->Factz,2));
	}

	
	

	//printf("cm %1.5e\t %1.5e %1.5e\n", vesicle->cm[0], vesicle->cm[1], vesicle->cm[2]);
	//printf("cm0 %1.5e\t %1.5e %1.5e\n", x0, y0, z0);
	/*
	for (int i=0; i<cstlist->n; i++){
		x0+=xlist[i];
		y0+=ylist[i];
		z0+=zlist[i];
	}
	x0/=cstlist->n;
	y0/=cstlist->n;
	z0/=cstlist->n;
	*/

	ts_double axx, axy, azx, ayy, ayz, azz;
	axx=0; axy=0; ayz=0; azx=0; ayy=0; azz=0;

	ts_double A[3][3];

	for (i=0; i<cstlist->n; i++){
		axx+=wlist[i]*(xlist[i]-vesicle->cm[0])*(xlist[i]-vesicle->cm[0]);
		ayy+=wlist[i]*(ylist[i]-vesicle->cm[1])*(ylist[i]-vesicle->cm[1]);
		azz+=wlist[i]*(zlist[i]-vesicle->cm[2])*(zlist[i]-vesicle->cm[2]);
		axy+=wlist[i]*(xlist[i]-vesicle->cm[0])*(ylist[i]-vesicle->cm[1]);
		ayz+=wlist[i]*(ylist[i]-vesicle->cm[1])*(zlist[i]-vesicle->cm[2]);
		azx+=wlist[i]*(zlist[i]-vesicle->cm[2])*(xlist[i]-vesicle->cm[0]);
	}
	
	A[0][0]=axx;
	A[1][1]=ayy;
	A[2][2]=azz;

	A[0][1]=axy;
	A[1][0]=axy;

	A[1][2]=ayz;
	A[2][1]=ayz;
	
	A[2][0]=azx;
	A[0][2]=azx;
	find_vesicle_dir(vesicle, A, dir);
	cluster_list_free(cstlist);
	free_vtx_cluster(vesicle);

}

void calculate_vertex_normal(ts_vesicle *vesicle){
		int ii, jj;
		
		ts_double xnorm, ynorm, znorm, norml;
		

		for (jj=0; jj<vesicle->vlist->n; jj++){
			if (vesicle->vlist->vtx[jj]->c>1e-15){
				xnorm=0.0;
				ynorm=0.0;
				znorm=0.0;
				for (ii=0;ii<vesicle->vlist->vtx[jj]->tristar_no; ii++){
					xnorm+=vesicle->vlist->vtx[jj]->tristar[ii]->xnorm;
					ynorm+=vesicle->vlist->vtx[jj]->tristar[ii]->ynorm;
					znorm+=vesicle->vlist->vtx[jj]->tristar[ii]->znorm;
				}

				norml=sqrt(xnorm*xnorm+ynorm*ynorm+znorm*znorm);
				xnorm/=norml;
				ynorm/=norml;
				znorm/=norml;
				vesicle->vlist->vtx[jj]->Fx=-xnorm;
				vesicle->vlist->vtx[jj]->Fy=-ynorm;
				vesicle->vlist->vtx[jj]->Fz=-znorm;
				//double FF=sqrt(pow(vesicle->vlist->vtx[jj]->Factx,2)+pow(vesicle->vlist->vtx[jj]->Facty,2)+pow(vesicle->vlist->vtx[jj]->Factz,2));
				//printf("%1.5e\n", FF );
				if (fabs(vesicle->vlist->vtx[jj]->Factx)<1e-15 && fabs(vesicle->vlist->vtx[jj]->Facty)<1e-15 && fabs(vesicle->vlist->vtx[jj]->Factz)<1e-15){
					vesicle->vlist->vtx[jj]->Factx=-vesicle->tape->F*xnorm;
					vesicle->vlist->vtx[jj]->Facty=-vesicle->tape->F*ynorm;
					vesicle->vlist->vtx[jj]->Factz=-vesicle->tape->F*znorm;
				}
			}
		}
}


void calculate_conc(ts_vesicle *vesicle){
	ts_double cs, conc;
	cs=vesicle->tape->cs;
	double c0=vesicle->tape->conc0;
	int i;
	for(i=0;i<vesicle->vlist->n;i++){
		//vtx->proj=vesicle->xnorm*(vtx->x-vesicle->cm[0])+vesicle->ynorm*(vtx->y-vesicle->cm[1])+vesicle->znorm*(vtx->z-vesicle->cm[2]);
		conc=c0*vesicle->tape->beta*(vesicle->proj_max-vesicle->proj_min)*vesicle->vmag*exp(-vesicle->tape->beta*vesicle->vmag*(vesicle->vlist->vtx[i]->proj-vesicle->proj_min)/vesicle->tape->D)/(1-exp(-vesicle->tape->beta*vesicle->vmag*(vesicle->proj_max-vesicle->proj_min)/vesicle->tape->D))/vesicle->tape->D;
		conc/=vesicle->volume;
		vesicle->vlist->vtx[i]->cx=conc;
	}

}

void calculate_proj_max_min(ts_vesicle *vesicle){
	
	double proj;
	int i;
	vesicle->proj_max=-1e7;
	vesicle->proj_min=1e7;

	double normal, xnorm, ynorm, znorm;


	
	xnorm=vesicle->xnorm;
	ynorm=vesicle->ynorm;
	znorm=vesicle->znorm;
	
	//normal=sqrt(vesicle->xnorm*vesicle->xnorm+vesicle->ynorm*vesicle->ynorm);

	//xnorm=vesicle->xnorm/normal;
	//ynorm=vesicle->ynorm/normal;

	for(i=0;i<vesicle->vlist->n;i++){
		proj=xnorm*(vesicle->vlist->vtx[i]->x-vesicle->cm[0])+ynorm*(vesicle->vlist->vtx[i]->y-vesicle->cm[1])+znorm*(vesicle->vlist->vtx[i]->z-vesicle->cm[2]);
		vesicle->vlist->vtx[i]->proj=proj;
		if (proj>vesicle->proj_max){
			vesicle->proj_max=proj;
		}
		if (proj<vesicle->proj_min){
			vesicle->proj_min=proj;
		}
	}
	//printf("proj %1.4e\t%1.4e\n", vesicle->proj_min, vesicle->proj_max);
}

void calculate_xback(ts_vesicle *vesicle){
	double proj;
	int i;
	//vesicle->proj_max=-1e7;
	vesicle->xback=1e7;

	for(i=0;i<vesicle->vlist->n;i++){
		if (vesicle->vlist->vtx[i]->x<vesicle->xback){
			vesicle->xback=vesicle->vlist->vtx[i]->x;
		}
	}
}




void UCSP_step(ts_vesicle *vesicle){
	ts_uint i, j;
	ts_double norml, proj, norml_old, xnorm_old, ynorm_old, znorm_old;
	ts_double cs, conc;
	cs=vesicle->tape->cs;
	double c0=vesicle->tape->conc0;
	double dir[3];
	ts_double F=vesicle->tape->F;
	int count=0;
	norml_old=0;
	xnorm_old=0;
	ynorm_old=0;
	znorm_old=0;
	FILE *fptr, *fptr1, *fptr2;
	
	
	

	
	while (1){
		calculate_proj_max_min(vesicle);
		for (j=0; j<vesicle->vlist->n; j++){
				proj=vesicle->vlist->vtx[j]->proj;
				//conc=c0*vesicle->tape->beta*vesicle->vmag*exp(-vesicle->tape->beta*vesicle->vmag*(proj-vesicle->proj_min)/vesicle->tape->D)/(1-exp(-vesicle->tape->beta*vesicle->vmag*(vesicle->proj_max-vesicle->proj_min)/vesicle->tape->D))/vesicle->tape->D;
				conc=c0*vesicle->tape->beta*(vesicle->proj_max-vesicle->proj_min)*vesicle->vmag*exp(-vesicle->tape->beta*vesicle->vmag*(proj-vesicle->proj_min)/vesicle->tape->D)/(1-exp(-vesicle->tape->beta*vesicle->vmag*(vesicle->proj_max-vesicle->proj_min)/vesicle->tape->D))/vesicle->tape->D;
				conc/=vesicle->volume;
				F=vesicle->vlist->vtx[j]->F_vtx*cs/(conc+cs);
				vesicle->vlist->vtx[j]->cx=conc;
				vesicle->vlist->vtx[j]->Factx=vesicle->vlist->vtx[j]->Fx*F;
				vesicle->vlist->vtx[j]->Facty=vesicle->vlist->vtx[j]->Fy*F;
				vesicle->vlist->vtx[j]->Factz=vesicle->vlist->vtx[j]->Fz*F;			
		}
		calculate_vesicle_dir(vesicle, dir);
		vesicle->xnorm=dir[0];
		vesicle->ynorm=dir[1];
		vesicle->znorm=dir[2];
		

		calculate_velocity_mag(vesicle, &vesicle->vmag);
		
		//printf("%1.5f\t%1.5f\t%1.5f\n", vesicle->xnorm, vesicle->ynorm, vesicle->znorm);

		double err;
		
		//printf("dir old : %1.4e\t%1.4e\t%1.4e\n", xnorm_old, ynorm_old, znorm_old);
		//printf("dir : %1.4e\t%1.4e\t%1.4e\n", vesicle->xnorm*vesicle->vmag, vesicle->ynorm*vesicle->vmag, vesicle->znorm*vesicle->vmag);
		err=(pow(vesicle->xnorm*vesicle->vmag-xnorm_old, 2)+pow(vesicle->ynorm*vesicle->vmag-ynorm_old, 2)+pow(vesicle->znorm*vesicle->vmag-znorm_old, 2))/(pow(vesicle->xnorm*vesicle->vmag, 2)+pow(vesicle->ynorm*vesicle->vmag, 2)+pow(vesicle->znorm*vesicle->vmag, 2));
		//printf("error -> %d\t%e\t%e\t%e\t%e\t%e\n", count, err, vesicle->xnorm, vesicle->ynorm, vesicle->znorm, vesicle->vmag);

		if (err<1e-3){
			//printf("converged in %d steps\n",count);
			break;
		}
		else{
			xnorm_old=vesicle->xnorm*vesicle->vmag;
			ynorm_old=vesicle->ynorm*vesicle->vmag;
			znorm_old=vesicle->znorm*vesicle->vmag;
		}
		count+=1;
		if (count>100000){
			printf("Converegence is too slow\n");
			exit(1);
			break;
		}
	}

	//printf("%1.5e \t %1.5e \t %1.5e \n", vesicle->xnorm, vesicle->ynorm, vesicle->znorm);

	fptr = fopen("./UCSP_step.d", "w");
	fprintf(fptr, "%1.5e \t %1.5e \t %1.5e \n", vesicle->xnorm, vesicle->ynorm, vesicle->znorm);
	fclose(fptr);	

	fptr1 = fopen("./conc.d", "w");
	for (j=0; j<vesicle->vlist->n; j++){
			//if (vesicle->vlist->vtx[j]->c!=0){
				proj=vesicle->vlist->vtx[j]->proj;
				//conc=c0*vesicle->tape->beta*vesicle->vmag*exp(-vesicle->tape->beta*vesicle->vmag*(proj-vesicle->proj_min)/vesicle->tape->D)/(1-exp(-vesicle->tape->beta*vesicle->vmag*(vesicle->proj_max-vesicle->proj_min)/vesicle->tape->D))/vesicle->tape->D;
				conc=c0*vesicle->tape->beta*(vesicle->proj_max-vesicle->proj_min)*vesicle->vmag*exp(-vesicle->tape->beta*vesicle->vmag*(proj-vesicle->proj_min)/vesicle->tape->D)/(1-exp(-vesicle->tape->beta*vesicle->vmag*(vesicle->proj_max-vesicle->proj_min)/vesicle->tape->D))/vesicle->tape->D;
				conc/=vesicle->volume;
				F= vesicle->vlist->vtx[j]->F_vtx*cs/(conc+cs);
				fprintf(fptr1, "%1.5e \t %1.5e \n", proj, conc);
				vesicle->vlist->vtx[j]->Factx=vesicle->vlist->vtx[j]->Fx*F;
				vesicle->vlist->vtx[j]->Facty=vesicle->vlist->vtx[j]->Fy*F;
				vesicle->vlist->vtx[j]->Factz=vesicle->vlist->vtx[j]->Fz*F;
			//}
	}
	fclose(fptr1);
	fptr2 = fopen("./speed.d", "w");
	fprintf(fptr2, "%1.5e %1.5e\n", vesicle->vmag, vesicle->volume);
	fclose(fptr2);

}
/*Shubhadeep*/



ts_bool run_simulation(ts_vesicle *vesicle, ts_uint mcsweeps, ts_uint inititer, ts_uint iterations, ts_uint start_iteration){
	ts_uint i, j, k, ii, jj; //,l,m;
	ts_double kc1=0,kc2=0,kc3=0,kc4=0;
	ts_double l1,l2,l3,vmsr,bfsr, vmsrt, bfsrt;
	ts_ulong epochtime;
	ts_double max_z;
	FILE *fd3=NULL;

	FILE *pFile2;

 	char filename[10000];
	//struct stat st;
	strcpy(filename,command_line_args.path);
	strcat(filename,"statistics.csv");
	
	//int result = stat(filename, &st);
	FILE *fd;
	if(start_iteration==0){
		fd=fopen(filename,"w");
		pFile2=fopen("global.txt", "w");
	}
	else{
		fd=fopen(filename,"a");
		pFile2=fopen("global.txt", "a");
	}
	if(fd==NULL){
		fatal("Cannot open statistics.csv file for writing",1);
	}
	if(start_iteration==0)
		fprintf(fd, "Epoch OuterLoop VertexMoveSucessRate BondFlipSuccessRate Volume Area lamdba1 lambda2 lambda3 Kc(2-9) Kc(6-9) Kc(2-end) Kc(3-6)\n");

/*	 if(vesicle->sphHarmonics!=NULL){
        strcpy(filename,command_line_args.path);
        strcat(filename,"ulm2.csv"); 
//	int result = stat(filename, &st);
	if(start_iteration==0)
		fd2=fopen(filename,"w");
	else
		fd2=fopen(filename,"a");
	if(fd2==NULL){
		fatal("Cannot open ulm2.csv file for writing",1);
	} 

	if(start_iteration==0) //file does not exist
		fprintf(fd2, "Timestep u_00^2 u_10^2 u_11^2 u_20^2 ...\n");	
	}
*/

/* RANDOM SEED SET BY CURRENT TIME */
	
	epochtime=get_epoch();			
	srand48(epochtime);
	if(vesicle->tape->allow_xy_plane_movement==0){
		centermass(vesicle, start_iteration-1, pFile2);
	}
	else{
		compute_cm(vesicle);
	}

	//printf("%1.4e\t %1.4e\t%1.4e", vesicle->cm[0], vesicle->cm[1], vesicle->cm[2]);
	vesicle->clist->c_cm[0]=vesicle->cm[0];
	vesicle->clist->c_cm[1]=vesicle->cm[1];
	vesicle->clist->c_cm[2]=vesicle->cm[2];

	cell_occupation(vesicle);
	vesicle_volume(vesicle); //needed for constant volume at this moment
	vesicle_area(vesicle); //needed for constant area at this moment
	if(V0<0.000001) 
		V0=vesicle->volume; 
	ts_fprintf(stdout,"Setting volume V0=%.17f\n",V0);
	if(A0<0.000001)
		A0=vesicle->area;
	ts_fprintf(stdout,"Setting area A0=%.17f\n",A0);
	epsvol=4.0*sqrt(2.0*M_PI)/pow(3.0,3.0/4.0)*V0/pow(vesicle->tlist->n,3.0/2.0);
//	printf("epsvol=%e\n",epsvol);
	epsarea=A0/(ts_double)vesicle->tlist->n;

	int ucsp_count=0;
	
	if(start_iteration<inititer) ts_fprintf(stdout, "Starting simulation (first %d x %d MC sweeps will not be recorded on disk)\n", inititer, mcsweeps);
	
	// calculate the direction of vesicle //
	/*Shubhadeep*/
	/*
	for (j=0; j<vesicle->vlist->n; j++){
		if (vesicle->vlist->vtx[j]->c!=0){
			vesicle->xnorm+=vesicle->vlist->vtx[j]->Factx;
			vesicle->ynorm+=vesicle->vlist->vtx[j]->Facty;
			vesicle->znorm+=vesicle->vlist->vtx[j]->Factz;
		}
	}
	ts_double norml=sqrt(vesicle->xnorm*vesicle->xnorm+vesicle->ynorm*vesicle->ynorm+vesicle->znorm*vesicle->znorm);
	vesicle->vmag=norml;
	vesicle->xnorm/=vesicle->vmag;
	vesicle->ynorm/=vesicle->vmag;
	vesicle->znorm/=vesicle->vmag;
	*/

	double dir[3];

	if (vesicle->tape->inhibition_switch && vesicle->tape->number_of_vertices_with_c0!=0 && vesicle->tape->F!=0){
		calculate_vesicle_dir(vesicle, dir);
	}

	vesicle->xnorm=dir[0];
	vesicle->ynorm=dir[1];
	vesicle->znorm=dir[2];
	calculate_vertex_normal(vesicle);
	calculate_velocity_mag(vesicle, &vesicle->vmag);
	vesicle->vmag/=10.;
	
	if (!vesicle->tape->inhibition_switch){
		for (i=0; i<vesicle->vlist->n; i++){
			vesicle->vlist->vtx[i]->cx=0;
		}
	}


	if (vesicle->tape->blow_switch){
		calculate_xback(vesicle);
	}
	int idx;
	if (vesicle->tape->blow_switch==0){
		for (idx=0; idx<vesicle->vlist->n; idx++){
			vesicle->vlist->vtx[idx]->Fbx=0;
			vesicle->vlist->vtx[idx]->Fby=0;
			vesicle->vlist->vtx[idx]->Fbz=0;
		}
	}
	//printf("%1.5e\t%1.5e\t%1.5e\n", vesicle->xnorm, vesicle->ynorm, vesicle->znorm );
	
	/*Shubhadeep*/
	for(i=start_iteration;i<inititer+iterations;i++){
		vmsr=0.0;
		bfsr=0.0;
		
	//plane confinement
	if(vesicle->tape->plane_confinement_switch){
		max_z=-1e10;
		for(k=0;k<vesicle->vlist->n;k++){
			if(vesicle->vlist->vtx[k]->z > max_z) max_z=vesicle->vlist->vtx[k]->z;
		}
		vesicle->confinement_plane.force_switch=0;
		if(max_z>=vesicle->tape->plane_d){
			ts_fprintf(stdout, "Max vertex out of bounds (z>=%e). Plane set to max_z = %e.\n",vesicle->tape->plane_d,max_z);
			vesicle->confinement_plane.z_max = max_z;
			vesicle->confinement_plane.force_switch=1;
		} else {
			vesicle->confinement_plane.z_max=vesicle->tape->plane_d;
		}

	    vesicle->confinement_plane.z_min=vesicle->tape->z_adhesion - 2*vesicle->tape->adhesion_radius;
		
		if(vesicle->confinement_plane.force_switch) ts_fprintf(stdout,"Squeezing with force %e.\n",vesicle->tape->plane_F);
	}
	//end plane confinement

//adhesion
	if(vesicle->tape->type_of_adhesion_model==3 || vesicle->tape->type_of_adhesion_model==4 || vesicle->tape->type_of_adhesion_model==5){	
		vesicle->adhesion_center = vesicle->tape->z_adhesion - vesicle->tape->adhesion_radius;
	}
//end of adhesion
	
/*    vesicle_volume(vesicle);
    fprintf(stderr,"Volume before TS=%1.16e\n", vesicle->volume); */

		
		for(j=0;j<mcsweeps;j++){
			double err=error_check(vesicle);

			if (vesicle->tape->F_noise_switch && (i*mcsweeps+j)%vesicle->tape->F_noise_interval==0){
				printf("Time step for random profiling: %d\n",(i*mcsweeps+j) );
				vesicle_random_force(vesicle);
			}

			

			//printf("before: %1.5e\t%1.5e\t%1.5e\n", vesicle->xnorm, vesicle->ynorm, vesicle->znorm);
			if (err>=0.1 && vesicle->tape->inhibition_switch && vesicle->tape->number_of_vertices_with_c0!=0 && vesicle->tape->F!=0){
				//printf("Calling UCSP --- %d, %1.5e\n", (j+1)*(i+1), err);
				UCSP_step(vesicle);
				ucsp_count+=1;
			}
			

			single_timestep(vesicle, &vmsrt, &bfsrt);

			/*
			printf("after steps\n" );
			for (idx=0; idx<vesicle->vlist->n; idx++){
				double Fmag=sqrt(pow(vesicle->vlist->vtx[idx]->Factx,2)+pow(vesicle->vlist->vtx[idx]->Facty,2)+pow(vesicle->vlist->vtx[idx]->Factz,2)) ;
				if (Fmag!=0){
					printf("Fx Fy Fz: %d\t%1.6e\t %1.6e\t%1.6e\t%1.6e\n",idx, vesicle->vlist->vtx[idx]->Factx, vesicle->vlist->vtx[idx]->Facty, vesicle->vlist->vtx[idx]->Factz, Fmag);
				}
			}
			*/
			vesicle_volume(vesicle); //calculates just volume. 
        	vesicle_area(vesicle); //calculates area.
			
			if (vesicle->tape->blow_switch){
				calculate_xback(vesicle);
			}
			
			compute_cm(vesicle);
			
			calculate_proj_max_min(vesicle);
			if (vesicle->tape->inhibition_switch){
				calculate_conc(vesicle);
			}
			
			//printf("after  %1.5e\t%1.5e\t%1.5e\n", vesicle->xnorm, vesicle->ynorm, vesicle->znorm);
			vmsr+=vmsrt;
			bfsr+=bfsrt;
		}
/*
    vesicle_volume(vesicle);
    fprintf(stderr,"Volume after TS=%1.16e\n", vesicle->volume); */
		vmsr/=(ts_double)mcsweeps;
		bfsr/=(ts_double)mcsweeps;
		if(vesicle->tape->allow_xy_plane_movement==0){
			centermass(vesicle, i, pFile2);
			vesicle->clist->c_cm[0]=vesicle->cm[0];
            vesicle->clist->c_cm[1]=vesicle->cm[1];
            vesicle->clist->c_cm[2]=vesicle->cm[2];
		}
		else{
			vesicle->clist->c_cm[0]=vesicle->cm[0];
            vesicle->clist->c_cm[1]=vesicle->cm[1];
            vesicle->clist->c_cm[2]=vesicle->cm[2];
		}
		cell_occupation(vesicle);
		
        dump_state(vesicle,i);
		
		if(vesicle->tape->constvolswitch==0){
			V0=vesicle->volume;
		}
		if(vesicle->tape->constareaswitch==0){
			A0=vesicle->area;
		}
		if(i>=inititer){
			write_vertex_xml_file(vesicle,i-inititer,NULL);
			write_master_xml_file(command_line_args.output_fullfilename);
			epochtime=get_epoch();			
			gyration_eigen(vesicle, &l1, &l2, &l3);
			//r0=getR0(vesicle);
/*            if(vesicle->sphHarmonics!=NULL){
			    preparationSh(vesicle,r0);
			    //calculateYlmi(vesicle);
			    calculateUlmComplex(vesicle);
			    storeUlmComplex2(vesicle);
			    saveAvgUlm2(vesicle);
                kc1=calculateKc(vesicle, 2,9);
                kc2=calculateKc(vesicle, 6,9);
                kc3=calculateKc(vesicle, 2,vesicle->sphHarmonics->l);
                kc4=calculateKc(vesicle, 3,6);

                strcpy(filename,command_line_args.path);
                strcat(filename,"state.dat");  
				fd1=fopen(filename,"w");
				fprintf(fd1,"%e %e\n",vesicle->volume, getR0(vesicle));
				for(k=0;k<vesicle->vlist->n;k++){
					fprintf(fd1,"%e %e %e %e %e\n",
						vesicle->vlist->vtx[k]->x,
						vesicle->vlist->vtx[k]->y,
						vesicle->vlist->vtx[k]->z,
						vesicle->vlist->vtx[k]->solAngle,
						vesicle->vlist->vtx[k]->relR
					);
				}
				fclose(fd1);
		
			fprintf(fd2,"%u ", i);
			for(l=0;l<vesicle->sphHarmonics->l;l++){
				for(m=l;m<2*l+1;m++){
					fprintf(fd2,"%e ", gsl_complex_abs2(vesicle->sphHarmonics->ulmComplex[l][m]) );
				}
			}
				fprintf(fd2,"\n");
	
		    	fflush(fd2);	


            }
*/

			fprintf(fd, "%lu %u %e %e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e\n",epochtime,i,vmsr,bfsr,vesicle->volume, vesicle->area,l1,l2,l3,kc1, kc2, kc3,kc4);

		    fflush(fd);	
		//	sprintf(filename,"timestep-%05d.pov",i-inititer);
		//	write_pov_file(vesicle,filename);
		} //end if(inititer....)
		fd3=fopen(".status","w"); //write status file when everything is written to disk.
		if(fd3==NULL){
			fatal("Cannot open .status file for writing",1);
		}
		fprintf(fd3,"%d",i);
		fclose(fd3);
		ts_fprintf(stdout,"Done %d out of %d iterations (x %d MC sweeps).\n",i+1,inititer+iterations,mcsweeps);

	}
	fclose(fd);
	fclose(pFile2);
//	if(fd2!=NULL) fclose(fd2);
	return TS_SUCCESS;
}

ts_bool single_timestep(ts_vesicle *vesicle,ts_double *vmsr, ts_double *bfsr){
//    vesicle_volume(vesicle);
//    fprintf(stderr,"Volume before TS=%1.16e\n", vesicle->volume);
    ts_bool retval;
    ts_double rnvec[3];
    ts_uint i,j, b;
    ts_uint vmsrcnt=0;

    for(i=0;i<vesicle->vlist->n;i++){
        rnvec[0]=drand48();
        rnvec[1]=drand48();
        rnvec[2]=drand48();
        retval=single_verticle_timestep(vesicle,vesicle->vlist->vtx[i],rnvec);
	if(retval==TS_SUCCESS) vmsrcnt++;        
    }


	ts_int bfsrcnt=0;
    for(i=0;i<3*vesicle->vlist->n;i++){
	b=rand() % vesicle->blist->n;
        //find a bond and return a pointer to a bond...
        //call single_bondflip_timestep...
        retval=single_bondflip_timestep(vesicle,vesicle->blist->bond[b],rnvec);
       // b++; retval=TS_FAIL;
	if(retval==TS_SUCCESS) bfsrcnt++;        
    }

	for(i=0;i<vesicle->poly_list->n;i++){
		for(j=0;j<vesicle->poly_list->poly[i]->vlist->n;j++){
			rnvec[0]=drand48();
			rnvec[1]=drand48();
			rnvec[2]=drand48();
			retval=single_poly_vertex_move(vesicle,vesicle->poly_list->poly[i],vesicle->poly_list->poly[i]->vlist->vtx[j],rnvec);	
		}
	}


	for(i=0;i<vesicle->filament_list->n;i++){
		for(j=0;j<vesicle->filament_list->poly[i]->vlist->n;j++){
			rnvec[0]=drand48();
			rnvec[1]=drand48();
			rnvec[2]=drand48();
			retval=single_filament_vertex_move(vesicle,vesicle->filament_list->poly[i],vesicle->filament_list->poly[i]->vlist->vtx[j],rnvec);	
		}
	}
 


//	printf("Bondflip success rate in one sweep: %d/%d=%e\n", cnt,3*vesicle->blist->n,(double)cnt/(double)vesicle->blist->n/3.0);
	*vmsr=(ts_double)vmsrcnt/(ts_double)vesicle->vlist->n;
	*bfsr=(ts_double)bfsrcnt/(ts_double)vesicle->vlist->n/3.0;
//    vesicle_volume(vesicle);
//    fprintf(stderr,"Volume after TS=%1.16e\n", vesicle->volume);
    return TS_SUCCESS;
}


    
