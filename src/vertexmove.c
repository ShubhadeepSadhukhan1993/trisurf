/* vim: set ts=4 sts=4 sw=4 noet : */
#include<stdlib.h>
#include<math.h>
#include "general.h"
#include "vertex.h"
#include "bond.h"
#include "triangle.h"
#include "vesicle.h"
#include "energy.h"
#include "timestep.h"
#include "cell.h"
//#include "io.h"
#include "io.h"
#include<stdio.h>
#include "vertexmove.h"
#include <string.h>
#include "constvol.h"

ts_bool single_verticle_timestep(ts_vesicle *vesicle,ts_vertex *vtx,ts_double *rn){
    ts_uint i;
    ts_double dist;
    ts_bool retval; 
    ts_uint cellidx; 
    ts_double delta_energy, delta_energy_cv,oenergy,dvol=0.0, darea=0.0, dstretchenergy=0.0;
    ts_double costheta,sintheta,phi,r;
	//This will hold all the information of vtx and its neighbours
	ts_vertex backupvtx[20], *constvol_vtx_moved=NULL, *constvol_vtx_backup=NULL;
	memcpy((void *)&backupvtx[0],(void *)vtx,sizeof(ts_vertex));

	//Some stupid tests for debugging cell occupation!
/*     	cellidx=vertex_self_avoidance(vesicle, vtx);
	if(vesicle->clist->cell[cellidx]==vtx->cell){
		fprintf(stderr,"Idx match!\n");
	} else {
		fprintf(stderr,"***** Idx don't match!\n");
		fatal("ENding.",1);
	}
*/

    	//temporarly moving the vertex
//	vtx->x=vtx->x+vesicle->stepsize*(2.0*rn[0]-1.0);
//    	vtx->y=vtx->y+vesicle->stepsize*(2.0*rn[1]-1.0);
//    	vtx->z=vtx->z+vesicle->stepsize*(2.0*rn[2]-1.0);

//random move in a sphere with radius stepsize:
	r=vesicle->stepsize*rn[0];
	phi=rn[1]*2*M_PI;
	costheta=2*rn[2]-1;
	sintheta=sqrt(1-pow(costheta,2));
	vtx->x=vtx->x+r*sintheta*cos(phi);
	vtx->y=vtx->y+r*sintheta*sin(phi);
	vtx->z=vtx->z+r*costheta;


//distance with neighbours check
    for(i=0;i<vtx->neigh_no;i++){
        dist=vtx_distance_sq(vtx,vtx->neigh[i]);
        if(dist<1.0 || dist>vesicle->dmax) {
		    vtx=memcpy((void *)vtx,(void *)&backupvtx[0],sizeof(ts_vertex));
    		return TS_FAIL;
		}
    }

// Distance with grafted poly-vertex check:	
	if(vtx->grafted_poly!=NULL){
		dist=vtx_distance_sq(vtx,vtx->grafted_poly->vlist->vtx[0]);
        if(dist<1.0 || dist>vesicle->dmax) {
		vtx=memcpy((void *)vtx,(void *)&backupvtx[0],sizeof(ts_vertex));
		return TS_FAIL;
		}
	}

// TODO: Maybe faster if checks only nucleus-neighboring cells
// Nucleus penetration check:
// #define SQ(x) x*x
if(vesicle->R_nucleus>0.0){
	if ((vtx->x-vesicle->nucleus_center[0])*(vtx->x-vesicle->nucleus_center[0])+ (vtx->y-vesicle->nucleus_center[1])*(vtx->y-vesicle->nucleus_center[1]) + (vtx->z-vesicle->nucleus_center[2])*(vtx->z-vesicle->nucleus_center[2]) < vesicle->R_nucleus){
		vtx=memcpy((void *)vtx,(void *)&backupvtx[0],sizeof(ts_vertex));
		return TS_FAIL;
	}
} else if(vesicle->R_nucleusX>0.0){
//	fprintf(stderr,"DEBUG, (Rx, Ry,Rz)^2=(%f,%f,%f)\n",vesicle->R_nucleusX, vesicle->R_nucleusY, vesicle->R_nucleusZ);
//	if (SQ(vtx->x-vesicle->nucleus_center[0])/vesicle->R_nucleusX + SQ(vtx->y-vesicle->nucleus_center[1])/vesicle->R_nucleusY + SQ(vtx->z-vesicle->nucleus_center[2])/vesicle->R_nucleusZ < 1.0){
	if ((vtx->x-vesicle->nucleus_center[0])*(vtx->x-vesicle->nucleus_center[0])/vesicle->R_nucleusX + (vtx->y-vesicle->nucleus_center[1])*(vtx->y-vesicle->nucleus_center[1])/vesicle->R_nucleusY + (vtx->z-vesicle->nucleus_center[2])*(vtx->z-vesicle->nucleus_center[2])/vesicle->R_nucleusZ < 1.0){
//	if (SQ(vtx->x)/vesicle->R_nucleusX + SQ(vtx->y)/vesicle->R_nucleusY + SQ(vtx->z)/vesicle->R_nucleusZ < 1.0){
		vtx=memcpy((void *)vtx,(void *)&backupvtx[0],sizeof(ts_vertex));
		return TS_FAIL;
	}

}

	// plane confinement check whether the new position of vertex will be out of bounds
	if(vesicle->tape->plane_confinement_switch){
		if(vtx->z>vesicle->confinement_plane.z_max || vtx->z<vesicle->confinement_plane.z_min){
		vtx=memcpy((void *)vtx,(void *)&backupvtx[0],sizeof(ts_vertex));
		return TS_FAIL;
		}

	}

	// adhesion check whether the new position of vertex will be out of bounds
	if(vesicle->tape->adhesion_switch){
		if (vesicle->tape->type_of_adhesion_model==1 || vesicle->tape->type_of_adhesion_model==2 || vesicle->tape->type_of_adhesion_model==6 || vesicle->tape->type_of_adhesion_model==7 || vesicle->tape->type_of_adhesion_model==8 || vesicle->tape->type_of_adhesion_model==9){
			if(vtx->z<vesicle->tape->z_adhesion){
			vtx=memcpy((void *)vtx,(void *)&backupvtx[0],sizeof(ts_vertex));
			return TS_FAIL;
			}
		}
		if (vesicle->tape->type_of_adhesion_model==3){
    		if((pow(vesicle->adhesion_center - vtx->z,2) + pow(vtx->x,2) + pow(vtx->y,2)) < pow(vesicle->tape->adhesion_radius,2)){
			vtx=memcpy((void *)vtx,(void *)&backupvtx[0],sizeof(ts_vertex));
			return TS_FAIL;
			}
        }
		if (vesicle->tape->type_of_adhesion_model==4){
    		if((pow(vesicle->adhesion_center - vtx->z,2) + pow(vtx->x,2)) < pow(vesicle->tape->adhesion_radius,2)){
			vtx=memcpy((void *)vtx,(void *)&backupvtx[0],sizeof(ts_vertex));
			return TS_FAIL;
			}
		}

		if (vesicle->tape->type_of_adhesion_model==5){
    		if((pow(vesicle->adhesion_center - vtx->z,2) + pow(vtx->y,2)) < pow(vesicle->tape->adhesion_radius,2)){
			vtx=memcpy((void *)vtx,(void *)&backupvtx[0],sizeof(ts_vertex));
			return TS_FAIL;
			}
		}
	}

//#undef SQ
//self avoidance check with distant vertices
	cellidx=vertex_self_avoidance(vesicle, vtx);
	//check occupation number
	retval=cell_occupation_number_and_internal_proximity(vesicle->clist,cellidx,vtx);

    if(retval==TS_FAIL){
		vtx=memcpy((void *)vtx,(void *)&backupvtx[0],sizeof(ts_vertex));
        return TS_FAIL;
    } 
   
 
//if all the tests are successful, then energy for vtx and neighbours is calculated
	for(i=0;i<vtx->neigh_no;i++){
	memcpy((void *)&backupvtx[i+1],(void *)vtx->neigh[i],sizeof(ts_vertex));
	}

	if(vesicle->pswitch == 1 || vesicle->tape->constvolswitch>0){
		for(i=0;i<vtx->tristar_no;i++) dvol-=vtx->tristar[i]->volume;
	}

    if(vesicle->tape->constareaswitch==2){
		for(i=0;i<vtx->tristar_no;i++) darea-=vtx->tristar[i]->area;
    
    }
	//stretching energy 1 of 3
	if(vesicle->tape->stretchswitch==1){
		for(i=0;i<vtx->tristar_no;i++) dstretchenergy-=vtx->tristar[i]->energy;
	}
    delta_energy=0;


//    vesicle_volume(vesicle);
//    fprintf(stderr,"Volume in the beginning=%1.16e\n", vesicle->volume);

    //update the normals of triangles that share bead i.
    for(i=0;i<vtx->tristar_no;i++) triangle_normal_vector(vtx->tristar[i]);
	oenergy=vtx->energy;
    energy_vertex(vtx);
    delta_energy=vtx->xk*(vtx->energy - oenergy);
    //the same is done for neighbouring vertices
    for(i=0;i<vtx->neigh_no;i++){
        oenergy=vtx->neigh[i]->energy;
        energy_vertex(vtx->neigh[i]);
        delta_energy+=vtx->neigh[i]->xk*(vtx->neigh[i]->energy-oenergy);
    }

	if(vesicle->pswitch == 1 || vesicle->tape->constvolswitch >0){
		for(i=0;i<vtx->tristar_no;i++) dvol+=vtx->tristar[i]->volume;
        if(vesicle->pswitch==1) delta_energy-=vesicle->pressure*dvol;
	};

    if(vesicle->tape->constareaswitch==2){
        /* check whether the darea is gt epsarea */
		for(i=0;i<vtx->tristar_no;i++) darea+=vtx->tristar[i]->area;
        if(fabs(vesicle->area+darea-A0)>epsarea){
	        //restore old state.
 			vtx=memcpy((void *)vtx,(void *)&backupvtx[0],sizeof(ts_vertex));
	        	for(i=0;i<vtx->neigh_no;i++){
		        	vtx->neigh[i]=memcpy((void *)vtx->neigh[i],(void *)&backupvtx[i+1],sizeof(ts_vertex));
	        	}
            		for(i=0;i<vtx->tristar_no;i++) triangle_normal_vector(vtx->tristar[i]); 
            		//fprintf(stderr,"fajlam!\n");
            		return TS_FAIL;
		}


    }

	if(vesicle->tape->constvolswitch==2){
		/*check whether the dvol is gt than epsvol */
			//fprintf(stderr,"DVOL=%1.16e\n",dvol);
		if(fabs(vesicle->volume+dvol-V0)>epsvol){
			//restore old state.
 			vtx=memcpy((void *)vtx,(void *)&backupvtx[0],sizeof(ts_vertex));
	        	for(i=0;i<vtx->neigh_no;i++){
		        	vtx->neigh[i]=memcpy((void *)vtx->neigh[i],(void *)&backupvtx[i+1],sizeof(ts_vertex));
	        	}
            		for(i=0;i<vtx->tristar_no;i++) triangle_normal_vector(vtx->tristar[i]); 
            		//fprintf(stderr,"fajlam!\n");
            		return TS_FAIL;
		}

	} else
//    vesicle_volume(vesicle);
//    fprintf(stderr,"Volume before=%1.16e\n", vesicle->volume);
   if(vesicle->tape->constvolswitch == 1){
        retval=constvolume(vesicle, vtx, -dvol, &delta_energy_cv, &constvol_vtx_moved,&constvol_vtx_backup);
        if(retval==TS_FAIL){ // if we couldn't move the vertex to assure constant volume
            vtx=memcpy((void *)vtx,(void *)&backupvtx[0],sizeof(ts_vertex));
	        for(i=0;i<vtx->neigh_no;i++){
		        vtx->neigh[i]=memcpy((void *)vtx->neigh[i],(void *)&backupvtx[i+1],sizeof(ts_vertex));
	        }
            for(i=0;i<vtx->tristar_no;i++) triangle_normal_vector(vtx->tristar[i]); 
 //           fprintf(stderr,"fajlam!\n");
            return TS_FAIL;
        }
//    vesicle_volume(vesicle);
//    fprintf(stderr,"Volume after=%1.16e\n", vesicle->volume);
//    fprintf(stderr,"Volume after-dvol=%1.16e\n", vesicle->volume-dvol);
//    fprintf(stderr,"Denergy before=%e\n",delta_energy);
    
    delta_energy+=delta_energy_cv;
//    fprintf(stderr,"Denergy after=%e\n",delta_energy);
    }
/* Vertices with spontaneous curvature may have spontaneous force perpendicular to the surface of the vesicle. additional delta energy is calculated in this function */

    /* Shubhadeep */
    if (vesicle->tape->shear_switch){
		delta_energy+=shear_force_energy(vesicle,vtx,backupvtx);
	}

	/* Shubhadeep */

	if (vesicle->tape->force_balance_along_z_axis==0){
		delta_energy+=direct_force_energy(vesicle,vtx,backupvtx);
	}
	else if (vesicle->tape->force_balance_along_z_axis==1){
		delta_energy+=direct_force_energy_with_Fz_balance(vesicle,vtx,backupvtx);
	}
	
	

	//stretching energy 2 of 3
	if(vesicle->tape->stretchswitch==1){
		for(i=0;i<vtx->tristar_no;i++){ 
			stretchenergy(vesicle, vtx->tristar[i]);
			dstretchenergy+=vtx->tristar[i]->energy;
			}
	}

	delta_energy+=dstretchenergy;	
		
/* No poly-bond energy for now!
	if(vtx->grafted_poly!=NULL){
		delta_energy+=
			(pow(sqrt(vtx_distance_sq(vtx, vtx->grafted_poly->vlist->vtx[0])-1),2)-
			pow(sqrt(vtx_distance_sq(&backupvtx[0], vtx->grafted_poly->vlist->vtx[0])-1),2)) *vtx->grafted_poly->k;
	}
*/

// plane confinement energy due to compressing force
	if(vesicle->tape->plane_confinement_switch){
		if(vesicle->confinement_plane.force_switch){
			//substract old energy
			if(abs(vesicle->tape->plane_d-vesicle->confinement_plane.z_max)>1e-10) {
				delta_energy-=vesicle->tape->plane_F / pow(vesicle->confinement_plane.z_max-backupvtx[0].z,2);
				delta_energy+=vesicle->tape->plane_F / pow(vesicle->confinement_plane.z_max-vtx->z,2);
			}
		}
	}

// change in energy due to adhesion
if(vesicle->tape->adhesion_switch){
//1 for step potential
	if(vesicle->tape->type_of_adhesion_model==1){

		if(vtx->z-vesicle->tape->z_adhesion<vesicle->tape->adhesion_cuttoff){
				delta_energy-=vesicle->tape->adhesion_strength;
		}
		if(backupvtx[0].z-vesicle->tape->z_adhesion<vesicle->tape->adhesion_cuttoff){
				delta_energy+=vesicle->tape->adhesion_strength;
		}
	}
//2 for parabolic potential
	else if(vesicle->tape->type_of_adhesion_model==2){

		if((vtx->z-vesicle->tape->z_adhesion)<=vesicle->tape->adhesion_cuttoff && (backupvtx[0].z-vesicle->tape->z_adhesion)<=vesicle->tape->adhesion_cuttoff){
				delta_energy-=(vesicle->tape->adhesion_strength/pow(vesicle->tape->adhesion_cuttoff,2))*(vtx->z - backupvtx[0].z)*(vtx->z + backupvtx[0].z - 2*vesicle->tape->adhesion_cuttoff);
		}
		else if((vtx->z-vesicle->tape->z_adhesion)<=vesicle->tape->adhesion_cuttoff && (backupvtx[0].z-vesicle->tape->z_adhesion)>vesicle->tape->adhesion_cuttoff){
				delta_energy-=(vesicle->tape->adhesion_strength/pow(vesicle->tape->adhesion_cuttoff,2))*pow(vtx->z - vesicle->tape->adhesion_cuttoff,2);
		}
		else if((vtx->z-vesicle->tape->z_adhesion)>vesicle->tape->adhesion_cuttoff && (backupvtx[0].z-vesicle->tape->z_adhesion)<=vesicle->tape->adhesion_cuttoff){
				delta_energy+=(vesicle->tape->adhesion_strength/pow(vesicle->tape->adhesion_cuttoff,2))*pow(backupvtx[0].z - vesicle->tape->adhesion_cuttoff,2);
		}
	}
//3 for sphrerical adhesion substrate with constant potential
	else if(vesicle->tape->type_of_adhesion_model==3){
//		if(pow(vesicle->adhesion_center - vtx->z,2) + pow(vtx->x,2) + pow(vtx->y,2) < pow(vesicle->tape->adhesion_cuttoff + vesicle->tape->adhesion_radius,2)){
		if(pow(pow(vesicle->adhesion_center - vtx->z,2) + pow(vtx->x,2) + pow(vtx->y,2),0.5) - vesicle->tape->adhesion_radius < vesicle->tape->adhesion_cuttoff){
			delta_energy-=vesicle->tape->adhesion_strength;
		}
//		if(pow(vesicle->adhesion_center - backupvtx[0].z,2) + pow(backupvtx[0].x,2) + pow(backupvtx[0].y,2) < pow(vesicle->tape->adhesion_cuttoff + vesicle->tape->adhesion_radius,2)){
		if(pow(pow(vesicle->adhesion_center - backupvtx[0].z,2) + pow(backupvtx[0].x,2) + pow(backupvtx[0].y,2),0.5) - vesicle->tape->adhesion_radius < vesicle->tape->adhesion_cuttoff){
			delta_energy+=vesicle->tape->adhesion_strength;
		}
	}
//4 for cylindrical adhesive substrate with constant potential along y-axis
	else if(vesicle->tape->type_of_adhesion_model==4){
		if(pow(pow(vesicle->adhesion_center - vtx->z,2) + pow(vtx->x,2),0.5) - vesicle->tape->adhesion_radius < vesicle->tape->adhesion_cuttoff){
			delta_energy-=vesicle->tape->adhesion_strength;
		}
		if(pow(pow(vesicle->adhesion_center - backupvtx[0].z,2) + pow(backupvtx[0].x,2),0.5) - vesicle->tape->adhesion_radius < vesicle->tape->adhesion_cuttoff){
			delta_energy+=vesicle->tape->adhesion_strength;
		}
	}

//5 for cylindrical adhesive substrate with constant potential along x-axis
	else if(vesicle->tape->type_of_adhesion_model==5){
		if(pow(pow(vesicle->adhesion_center - vtx->z,2) + pow(vtx->y,2),0.5) - vesicle->tape->adhesion_radius < vesicle->tape->adhesion_cuttoff){
			delta_energy-=vesicle->tape->adhesion_strength;
		}
		if(pow(pow(vesicle->adhesion_center - backupvtx[0].z,2) + pow(backupvtx[0].y,2),0.5) - vesicle->tape->adhesion_radius < vesicle->tape->adhesion_cuttoff){
			delta_energy+=vesicle->tape->adhesion_strength;
		}
	}
//Shubhadeep
	else if(vesicle->tape->type_of_adhesion_model==6){
		ts_double x_max, x_min;
		ts_int i;
		x_max=-1e308;
		x_min=1e308;
		for (i=0; i<vesicle->vlist->n; i++){
			if (vesicle->vlist->vtx[i]->x > x_max){
				x_max=vesicle->vlist->vtx[i]->x;
			}
			if (vesicle->vlist->vtx[i]->x < x_min){
				x_min=vesicle->vlist->vtx[i]->x;
			}
		}
		//printf("xmax %1.2f\n", x_max);
		//printf("xmin %1.2f\n", x_min);
		//x_max*=1.5;
		//x_max*=1.5;
		ts_double dE_ad;
		dE_ad=vesicle->tape->adhesion_strength*(vtx->x-x_min)/(x_max-x_min); //vesicle->tape->adhesion_strength+2*
		//if (dE_ad>=vesicle->tape->adhesion_strength)
		//	printf("%1.5f\n", dE_ad);

		if(vtx->z-vesicle->tape->z_adhesion<vesicle->tape->adhesion_cuttoff){
				delta_energy-=dE_ad;
		}
		if(backupvtx[0].z-vesicle->tape->z_adhesion<vesicle->tape->adhesion_cuttoff){
				delta_energy+=dE_ad;
		}

	}
   //Shubhadeep
	else if(vesicle->tape->type_of_adhesion_model==7){
		//printf("Acting 7\n");
		if(vtx->z-vesicle->tape->z_adhesion<vesicle->tape->adhesion_cuttoff && vtx->c!=0){
			//printf("Acting 7 inside\n");
				delta_energy-=vesicle->tape->adhesion_strength;
		}
		if(backupvtx[0].z-vesicle->tape->z_adhesion<vesicle->tape->adhesion_cuttoff && vtx->c!=0){
			//printf("Acting 7 outside\n");
				delta_energy+=vesicle->tape->adhesion_strength;
		}

		//if(vtx->z-vesicle->tape->z_adhesion<vesicle->tape->adhesion_cuttoff && vtx->c!=0 && backupvtx[0].z-vesicle->tape->z_adhesion>vesicle->tape->adhesion_cuttoff){
		//	printf("done\n");
		//}
		
	}
	else if(vesicle->tape->type_of_adhesion_model==8){
		//printf("Acting 8 %1.3f\n",vesicle->tape->lamda);
		
		/*
		if(vtx->z-vesicle->tape->z_adhesion<vesicle->tape->adhesion_cuttoff){
			if ((int)(floor((vtx->x-vesicle->tape->lamda/2.)/vesicle->tape->lamda))%2!=0){
				delta_energy-=vesicle->tape->adhesion_strength;
			}
		}
		if(backupvtx[0].z-vesicle->tape->z_adhesion<vesicle->tape->adhesion_cuttoff){
			if ((int)(floor((backupvtx[0].x-vesicle->tape->lamda/2.)/vesicle->tape->lamda))%2!=0){
				delta_energy+=vesicle->tape->adhesion_strength;
			}
		}
		
		*/
		double xx;
		if(vtx->z-vesicle->tape->z_adhesion<vesicle->tape->adhesion_cuttoff){
			xx=(vtx->x+vesicle->tape->lamda1/2.)/(vesicle->tape->lamda1+vesicle->tape->lamda2);
			if (xx-floor(xx)<vesicle->tape->lamda1/(vesicle->tape->lamda1+vesicle->tape->lamda2)){
				delta_energy-=vesicle->tape->adhesion_strength;
			}
		}
		if(backupvtx[0].z-vesicle->tape->z_adhesion<vesicle->tape->adhesion_cuttoff){
			xx=(backupvtx[0].x+vesicle->tape->lamda1/2.)/(vesicle->tape->lamda1+vesicle->tape->lamda2);
			if (xx-floor(xx)<vesicle->tape->lamda1/(vesicle->tape->lamda1+vesicle->tape->lamda2)){
				delta_energy+=vesicle->tape->adhesion_strength;
			}
		}

		//if(vtx->z-vesicle->tape->z_adhesion<vesicle->tape->adhesion_cuttoff && vtx->c!=0 && backupvtx[0].z-vesicle->tape->z_adhesion>vesicle->tape->adhesion_cuttoff){
		//	printf("done\n");
		//}
		
	}
	else if(vesicle->tape->type_of_adhesion_model==9){
		ts_double x_max, x_min;
		ts_int i;
		x_max=-1e308;
		x_min=1e308;
		for (i=0; i<vesicle->vlist->n; i++){
			if (vesicle->vlist->vtx[i]->x > x_max){
				x_max=vesicle->vlist->vtx[i]->x;
			}
			if (vesicle->vlist->vtx[i]->x < x_min){
				x_min=vesicle->vlist->vtx[i]->x;
			}
		}
		//printf("xmax %1.2f\n", x_max);
		//printf("xmin %1.2f\n", x_min);
		//x_max*=1.5;
		//x_max*=1.5;
		ts_double dE_ad;
		dE_ad=2*vesicle->tape->adhesion_strength*fabs(vtx->x-x_min/2.-x_max/2.)/(x_max-x_min); //vesicle->tape->adhesion_strength+2*
		//if (dE_ad>=vesicle->tape->adhesion_strength)
		//	printf("%1.5f\n", dE_ad);

		if(vtx->z-vesicle->tape->z_adhesion<vesicle->tape->adhesion_cuttoff){
				delta_energy-=dE_ad;
		}
		if(backupvtx[0].z-vesicle->tape->z_adhesion<vesicle->tape->adhesion_cuttoff){
				delta_energy+=dE_ad;
		}
	}
//Shubhadeep
	
}

//   fprintf(stderr, "DE=%f\n",delta_energy);
    //MONTE CARLOOOOOOOO
//	if(vtx->c!=0.0) printf("DE=%f\n",delta_energy);
    if(delta_energy>=0){
#ifdef TS_DOUBLE_DOUBLE
        if(exp(-delta_energy)< drand48())
#endif
#ifdef TS_DOUBLE_FLOAT
        if(expf(-delta_energy)< (ts_float)drand48())
#endif
#ifdef TS_DOUBLE_LONGDOUBLE
        if(expl(-delta_energy)< (ts_ldouble)drand48())
#endif
    {
    //not accepted, reverting changes
  //  fprintf(stderr,"MC failed\n");
	vtx=memcpy((void *)vtx,(void *)&backupvtx[0],sizeof(ts_vertex));
	for(i=0;i<vtx->neigh_no;i++){
		vtx->neigh[i]=memcpy((void *)vtx->neigh[i],(void *)&backupvtx[i+1],sizeof(ts_vertex));
	}
	
    //update the normals of triangles that share bead i.
   for(i=0;i<vtx->tristar_no;i++) triangle_normal_vector(vtx->tristar[i]);
	//stretching energy 3 of 3
	if(vesicle->tape->stretchswitch==1){
		for(i=0;i<vtx->tristar_no;i++){ 
			stretchenergy(vesicle,vtx->tristar[i]);
			}
	}

//    fprintf(stderr, "before vtx(x,y,z)=%e,%e,%e\n",constvol_vtx_moved->x, constvol_vtx_moved->y, constvol_vtx_moved->z);
    if(vesicle->tape->constvolswitch == 1){
        constvolumerestore(constvol_vtx_moved,constvol_vtx_backup);
    }
//    fprintf(stderr, "after vtx(x,y,z)=%e,%e,%e\n",constvol_vtx_moved->x, constvol_vtx_moved->y, constvol_vtx_moved->z);
//    vesicle_volume(vesicle);
//    fprintf(stderr,"Volume after fail=%1.16e\n", vesicle->volume);
    return TS_FAIL; 
    }
}
	//accepted	
 //   fprintf(stderr,"MC accepted\n");
//	oldcellidx=vertex_self_avoidance(vesicle, &backupvtx[0]);
	if(vtx->cell!=vesicle->clist->cell[cellidx]){
		retval=cell_add_vertex(vesicle->clist->cell[cellidx],vtx);
//		if(retval==TS_SUCCESS) cell_remove_vertex(vesicle->clist->cell[oldcellidx],vtx);
		if(retval==TS_SUCCESS) cell_remove_vertex(backupvtx[0].cell,vtx);
		
	}

    if(vesicle->tape->constvolswitch == 2){
	vesicle->volume+=dvol;
    } else
    if(vesicle->tape->constvolswitch == 1){
        constvolumeaccept(vesicle,constvol_vtx_moved,constvol_vtx_backup);
    }

    if(vesicle->tape->constareaswitch==2){
        vesicle->area+=darea;
    }
//	if(oldcellidx);
    //END MONTE CARLOOOOOOO
//    vesicle_volume(vesicle);
//    fprintf(stderr,"Volume after success=%1.16e\n", vesicle->volume);
    return TS_SUCCESS;
}


ts_bool single_poly_vertex_move(ts_vesicle *vesicle,ts_poly *poly,ts_vertex *vtx,ts_double *rn){
	ts_uint i;
	ts_bool retval; 
	ts_uint cellidx; 
//	ts_double delta_energy;
	ts_double costheta,sintheta,phi,r;
	ts_double dist;
	//This will hold all the information of vtx and its neighbours
	ts_vertex backupvtx;
//	ts_bond backupbond[2];
	memcpy((void *)&backupvtx,(void *)vtx,sizeof(ts_vertex));

	//random move in a sphere with radius stepsize:
	r=vesicle->stepsize*rn[0];
	phi=rn[1]*2*M_PI;
	costheta=2*rn[2]-1;
	sintheta=sqrt(1-pow(costheta,2));
	vtx->x=vtx->x+r*sintheta*cos(phi);
	vtx->y=vtx->y+r*sintheta*sin(phi);
	vtx->z=vtx->z+r*costheta;


	//distance with neighbours check
	for(i=0;i<vtx->neigh_no;i++){
		dist=vtx_distance_sq(vtx,vtx->neigh[i]);
		if(dist<1.0 || dist>vesicle->dmax) {
			vtx=memcpy((void *)vtx,(void *)&backupvtx,sizeof(ts_vertex));
			return TS_FAIL;
		}
	}

// Distance with grafted vesicle-vertex check:	
	if(vtx==poly->vlist->vtx[0]){
		dist=vtx_distance_sq(vtx,poly->grafted_vtx);
        if(dist<1.0 || dist>vesicle->dmax) {
		vtx=memcpy((void *)vtx,(void *)&backupvtx,sizeof(ts_vertex));
		return TS_FAIL;
		}
	}


	//self avoidance check with distant vertices
	cellidx=vertex_self_avoidance(vesicle, vtx);
	//check occupation number
	retval=cell_occupation_number_and_internal_proximity(vesicle->clist,cellidx,vtx);
	
	if(retval==TS_FAIL){
		vtx=memcpy((void *)vtx,(void *)&backupvtx,sizeof(ts_vertex));
        return TS_FAIL;
	} 


	//if all the tests are successful, then energy for vtx and neighbours is calculated
/* Energy ignored for now!
	delta_energy=0;
	for(i=0;i<vtx->bond_no;i++){
		memcpy((void *)&backupbond[i],(void *)vtx->bond[i],sizeof(ts_bond));

		vtx->bond[i]->bond_length=sqrt(vtx_distance_sq(vtx->bond[i]->vtx1,vtx->bond[i]->vtx2));
		bond_energy(vtx->bond[i],poly);
		delta_energy+= vtx->bond[i]->energy - backupbond[i].energy;
	}

	if(vtx==poly->vlist->vtx[0]){
		delta_energy+=
			(pow(sqrt(vtx_distance_sq(vtx, poly->grafted_vtx)-1),2)-
			pow(sqrt(vtx_distance_sq(&backupvtx, poly->grafted_vtx)-1),2)) *poly->k;
		
	}


	if(delta_energy>=0){
#ifdef TS_DOUBLE_DOUBLE
        if(exp(-delta_energy)< drand48() )
#endif
#ifdef TS_DOUBLE_FLOAT
        if(expf(-delta_energy)< (ts_float)drand48())
#endif
#ifdef TS_DOUBLE_LONGDOUBLE
        if(expl(-delta_energy)< (ts_ldouble)drand48())
#endif
    	{
	//not accepted, reverting changes
	vtx=memcpy((void *)vtx,(void *)&backupvtx,sizeof(ts_vertex));
	for(i=0;i<vtx->bond_no;i++){
	vtx->bond[i]=memcpy((void *)vtx->bond[i],(void *)&backupbond[i],sizeof(ts_bond));
	}

    return TS_FAIL; 
	}
	}
*/
		
//	oldcellidx=vertex_self_avoidance(vesicle, &backupvtx[0]);
	if(vtx->cell!=vesicle->clist->cell[cellidx]){
		retval=cell_add_vertex(vesicle->clist->cell[cellidx],vtx);
//		if(retval==TS_SUCCESS) cell_remove_vertex(vesicle->clist->cell[oldcellidx],vtx);
		if(retval==TS_SUCCESS) cell_remove_vertex(backupvtx.cell,vtx);	
	}
//	if(oldcellidx);
    //END MONTE CARLOOOOOOO
    return TS_SUCCESS;
}




ts_bool single_filament_vertex_move(ts_vesicle *vesicle,ts_poly *poly,ts_vertex *vtx,ts_double *rn){
	ts_uint i;
	ts_bool retval; 
	ts_uint cellidx; 
	ts_double delta_energy;
	ts_double costheta,sintheta,phi,r;
	ts_double dist[2];
	//This will hold all the information of vtx and its neighbours
	ts_vertex backupvtx,backupneigh[2];
	ts_bond backupbond[2];

	//backup vertex:		
	memcpy((void *)&backupvtx,(void *)vtx,sizeof(ts_vertex));

	//random move in a sphere with radius stepsize:
	r=vesicle->stepsize*rn[0];
	phi=rn[1]*2*M_PI;
	costheta=2*rn[2]-1;
	sintheta=sqrt(1-pow(costheta,2));
	vtx->x=vtx->x+r*sintheta*cos(phi);
	vtx->y=vtx->y+r*sintheta*sin(phi);
	vtx->z=vtx->z+r*costheta;


	//distance with neighbours check
	for(i=0;i<vtx->bond_no;i++){
		dist[i]=vtx_distance_sq(vtx->bond[i]->vtx1,vtx->bond[i]->vtx2);
		if(dist[i]<1.0 || dist[i]>vesicle->dmax) {
			vtx=memcpy((void *)vtx,(void *)&backupvtx,sizeof(ts_vertex));
			return TS_FAIL;
		}
	}

// TODO: Maybe faster if checks only nucleus-neighboring cells
// Nucleus penetration check:
	if (vtx->x*vtx->x + vtx->y*vtx->y + vtx->z*vtx->z < vesicle->R_nucleus){
		vtx=memcpy((void *)vtx,(void *)&backupvtx,sizeof(ts_vertex));
		return TS_FAIL;
	}


	//self avoidance check with distant vertices
	cellidx=vertex_self_avoidance(vesicle, vtx);
	//check occupation number
	retval=cell_occupation_number_and_internal_proximity(vesicle->clist,cellidx,vtx);
	if(retval==TS_FAIL){
		vtx=memcpy((void *)vtx,(void *)&backupvtx,sizeof(ts_vertex));
        return TS_FAIL;
	} 

	//backup bonds
	for(i=0;i<vtx->bond_no;i++){
		memcpy(&backupbond[i],vtx->bond[i], sizeof(ts_bond));
		vtx->bond[i]->bond_length=sqrt(dist[i]);
		bond_vector(vtx->bond[i]);
	}

	//backup neighboring vertices:
	for(i=0;i<vtx->neigh_no;i++){
		memcpy(&backupneigh[i],vtx->neigh[i], sizeof(ts_vertex));
	}
	
	//if all the tests are successful, then energy for vtx and neighbours is calculated
	delta_energy=0;
	
	if(vtx->bond_no == 2){
		vtx->energy = -(vtx->bond[0]->x*vtx->bond[1]->x + vtx->bond[0]->y*vtx->bond[1]->y + vtx->bond[0]->z*vtx->bond[1]->z)/vtx->bond[0]->bond_length/vtx->bond[1]->bond_length;
		delta_energy += vtx->energy - backupvtx.energy;
	}

	for(i=0;i<vtx->neigh_no;i++){
		if(vtx->neigh[i]->bond_no == 2){
			vtx->neigh[i]->energy = -(vtx->neigh[i]->bond[0]->x*vtx->neigh[i]->bond[1]->x + vtx->neigh[i]->bond[0]->y*vtx->neigh[i]->bond[1]->y + vtx->neigh[i]->bond[0]->z*vtx->neigh[i]->bond[1]->z)/vtx->neigh[i]->bond[0]->bond_length/vtx->neigh[i]->bond[1]->bond_length;
			delta_energy += vtx->neigh[i]->energy - backupneigh[i].energy;
		}
	}

	// poly->k is filament persistence length (in units l_min)
	delta_energy *= poly->k;

	if(delta_energy>=0){
#ifdef TS_DOUBLE_DOUBLE
        if(exp(-delta_energy)< drand48() )
#endif
#ifdef TS_DOUBLE_FLOAT
        if(expf(-delta_energy)< (ts_float)drand48())
#endif
#ifdef TS_DOUBLE_LONGDOUBLE
        if(expl(-delta_energy)< (ts_ldouble)drand48())
#endif
    	{
	//not accepted, reverting changes
	vtx=memcpy((void *)vtx,(void *)&backupvtx,sizeof(ts_vertex));
	for(i=0;i<vtx->neigh_no;i++){
		memcpy(vtx->neigh[i],&backupneigh[i],sizeof(ts_vertex));
	}
	for(i=0;i<vtx->bond_no;i++){
		vtx->bond[i]=memcpy((void *)vtx->bond[i],(void *)&backupbond[i],sizeof(ts_bond));
	}

    return TS_FAIL; 
	}
	}
	
	
//	oldcellidx=vertex_self_avoidance(vesicle, &backupvtx[0]);
	if(vtx->cell!=vesicle->clist->cell[cellidx]){
		retval=cell_add_vertex(vesicle->clist->cell[cellidx],vtx);
//		if(retval==TS_SUCCESS) cell_remove_vertex(vesicle->clist->cell[oldcellidx],vtx);
		if(retval==TS_SUCCESS) cell_remove_vertex(backupvtx.cell,vtx);	
	}
//	if(oldcellidx);
    //END MONTE CARLOOOOOOO
    return TS_SUCCESS;
}
