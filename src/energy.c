/* vim: set ts=4 sts=4 sw=4 noet : */
#include<stdlib.h>
#include "general.h"
#include "energy.h"
#include "vertex.h"
#include<math.h>
#include<stdio.h>



// Shubhadeep //
const gsl_rng_type * TT;
gsl_rng * rr_noise;


FILE *fforce;


ts_double random_force(ts_vesicle *vesicle){
	return gsl_ran_gaussian_tail (rr_noise, -vesicle->tape->F, vesicle->tape->F_noise_SD)+vesicle->tape->F;
} 


void vesicle_random_force(ts_vesicle *vesicle){
	int i;
	fforce = fopen("./random_force.d", "a");
	for (i=0; i<vesicle->vlist->n; i++){
		if (vesicle->vlist->vtx[i]->c!=0){ 
			vesicle->vlist->vtx[i]->F_vtx=random_force(vesicle);
			fprintf(fforce, "%1.7e\t", vesicle->vlist->vtx[i]->F_vtx);
		}
		else{
			vesicle->vlist->vtx[i]->F_vtx=0;
		}
	}
	fprintf(fforce,"\n");
	fclose(fforce);	
}
// Shubhadeep //

//gsl_rng_env_setup();
//TT = gsl_rng_default;

/** @brief Wrapper that calculates energy of every vertex in vesicle
 *  
 *  Function calculated energy of every vertex in vesicle. It can be used in
 *  initialization procedure or in recalculation of the energy after non-MCsweep *  operations. However, when random move of vertex or flip of random bond occur *  call to this function is not necessary nor recommended. 
 *  @param *vesicle is a pointer to vesicle.
 *  @returns TS_SUCCESS on success.
*/
ts_bool mean_curvature_and_energy(ts_vesicle *vesicle){

    ts_uint i;
    
    ts_vertex_list *vlist=vesicle->vlist;
    ts_vertex **vtx=vlist->vtx;

    for(i=0;i<vlist->n;i++){
        energy_vertex(vtx[i]);
        
    }

    return TS_SUCCESS;
}

/** @brief Calculate energy of a bond (in models where energy is bond related)
 *
 *  This function is experimental and currently only used in polymeres calculation (PEGs or polymeres inside the vesicle).
 *
 *  @param *bond is a pointer to a bond between two vertices in polymere
 *  @param *poly is a pointer to polymere in which we calculate te energy of the bond
 *  @returns TS_SUCCESS on successful calculation
*/
inline ts_bool bond_energy(ts_bond *bond,ts_poly *poly){
//TODO: This value to be changed and implemented in data structure:
	ts_double d_relaxed=1.0;
	bond->energy=poly->k*pow(bond->bond_length-d_relaxed,2);
	return TS_SUCCESS;
};

/** @brief Calculation of the bending energy of the vertex.
 *  
 *  Main function that calculates energy of the vertex \f$i\f$. Function returns \f$\frac{1}{2}(c_1+c_2-c)^2 s\f$, where \f$(c_1+c_2)/2\f$ is mean curvature,
 * \f$c/2\f$ is spontaneous curvature and \f$s\f$ is area per vertex \f$i\f$.
 *
 * Nearest neighbors (NN) must be ordered in counterclockwise direction for this function to work.
 *  Firstly NNs that form two neighboring triangles are found (\f$j_m\f$, \f$j_p\f$ and common \f$j\f$). Later, the scalar product of vectors \f$x_1=(\mathbf{i}-\mathbf{j_p})\cdot (\mathbf{i}-\mathbf{j_p})(\mathbf{i}-\mathbf{j_p})\f$, \f$x_2=(\mathbf{j}-\mathbf{j_p})\cdot  (\mathbf{j}-\mathbf{j_p})\f$  and \f$x_3=(\mathbf{j}-\mathbf{j_p})\cdot (\mathbf{i}-\mathbf{j_p})\f$  are calculated. From these three vectors the \f$c_{tp}=\frac{1}{\tan(\varphi_p)}\f$ is calculated, where \f$\varphi_p\f$ is the inner angle at vertex \f$j_p\f$. The procedure is repeated for \f$j_m\f$ instead of \f$j_p\f$ resulting in \f$c_{tn}\f$.
 *  
\begin{tikzpicture}{
\coordinate[label=below:$i$] (i) at (2,0);
\coordinate[label=left:$j_m$] (jm) at (0,3.7);
\coordinate[label=above:$j$] (j) at (2.5,6.4);
\coordinate[label=right:$j_p$] (jp) at (4,2.7);

\draw (i) -- (jm) -- (j) -- (jp) -- (i) -- (j);

\begin{scope}
\path[clip] (jm)--(i)--(j);
\draw (jm) circle (0.8);
\node[right] at (jm) {$\varphi_m$};
\end{scope}

\begin{scope}
\path[clip] (jp)--(i)--(j);
\draw (jp) circle (0.8);
\node[left] at (jp) {$\varphi_p$};
\end{scope}

%%vertices
\draw [fill=gray] (i) circle (0.1);
\draw [fill=white] (j) circle (0.1);
\draw [fill=white] (jp) circle (0.1);
\draw [fill=white] (jm) circle (0.1);
%\node[draw,circle,fill=white] at (i) {};
\end{tikzpicture}

 * The curvature is then calculated as \f$\mathbf{h}=\frac{1}{2}\Sigma_{k=0}^{\mathrm{neigh\_no}} c_{tp}^{(k)}+c_{tm}^{(k)} (\mathbf{j_k}-\mathbf{i})\f$, where \f$c_{tp}^{(k)}+c_{tm}^k=2\sigma^{(k)}\f$ (length in dual lattice?) and the previous equation can be written as \f$\mathbf{h}=\Sigma_{k=0}^{\mathrm{neigh\_no}}\sigma^{(k)}\cdot(\mathbf{j}-\mathbf{i})\f$ (See Kroll, p. 384, eq 70).
 *
 * From the curvature the enery is calculated by equation \f$E=\frac{1}{2}\mathbf{h}\cdot\mathbf{h}\f$.
 * @param *vtx is a pointer to vertex at which we want to calculate the energy
 * @returns TS_SUCCESS on successful calculation.
*/

inline ts_bool energy_vertex(ts_vertex *vtx){
    ts_uint jj;
    ts_uint jjp,jjm;
    ts_vertex *j,*jp, *jm;
    ts_triangle *jt;
    ts_double s=0.0,xh=0.0,yh=0.0,zh=0.0,txn=0.0,tyn=0.0,tzn=0.0;
    ts_double x1,x2,x3,ctp,ctm,tot,xlen;
    ts_double h,ht;
    for(jj=1; jj<=vtx->neigh_no; jj++){
        jjp=jj+1;
        if(jjp>vtx->neigh_no) jjp=1;
        jjm=jj-1;
        if(jjm<1) jjm=vtx->neigh_no;
        j=vtx->neigh[jj-1];
        jp=vtx->neigh[jjp-1];
        jm=vtx->neigh[jjm-1];
        jt=vtx->tristar[jj-1];
        x1=vtx_distance_sq(vtx,jp); //shouldn't be zero!
        x2=vtx_distance_sq(j,jp);   // shouldn't be zero!
        x3=(j->x-jp->x)*(vtx->x-jp->x)+
           (j->y-jp->y)*(vtx->y-jp->y)+
           (j->z-jp->z)*(vtx->z-jp->z);
        
#ifdef TS_DOUBLE_DOUBLE
        ctp=x3/sqrt(x1*x2-x3*x3);
#endif
#ifdef TS_DOUBLE_FLOAT
        ctp=x3/sqrtf(x1*x2-x3*x3);
#endif
#ifdef TS_DOUBLE_LONGDOUBLE
        ctp=x3/sqrtl(x1*x2-x3*x3);
#endif
        x1=vtx_distance_sq(vtx,jm);
        x2=vtx_distance_sq(j,jm);
        x3=(j->x-jm->x)*(vtx->x-jm->x)+
           (j->y-jm->y)*(vtx->y-jm->y)+
           (j->z-jm->z)*(vtx->z-jm->z);
#ifdef TS_DOUBLE_DOUBLE
        ctm=x3/sqrt(x1*x2-x3*x3);
#endif
#ifdef TS_DOUBLE_FLOAT
        ctm=x3/sqrtf(x1*x2-x3*x3);
#endif
#ifdef TS_DOUBLE_LONGDOUBLE
        ctm=x3/sqrtl(x1*x2-x3*x3);
#endif
        tot=ctp+ctm;
        tot=0.5*tot;

        xlen=vtx_distance_sq(j,vtx);
/*
#ifdef  TS_DOUBLE_DOUBLE 
        vtx->bond[jj-1]->bond_length=sqrt(xlen); 
#endif
#ifdef  TS_DOUBLE_FLOAT
        vtx->bond[jj-1]->bond_length=sqrtf(xlen); 
#endif
#ifdef  TS_DOUBLE_LONGDOUBLE 
        vtx->bond[jj-1]->bond_length=sqrtl(xlen); 
#endif

        vtx->bond[jj-1]->bond_length_dual=tot*vtx->bond[jj-1]->bond_length;
*/
        s+=tot*xlen;
        xh+=tot*(j->x - vtx->x);
        yh+=tot*(j->y - vtx->y);
        zh+=tot*(j->z - vtx->z);
        txn+=jt->xnorm;
        tyn+=jt->ynorm;
        tzn+=jt->znorm;
    }
    
    h=xh*xh+yh*yh+zh*zh;
    ht=txn*xh+tyn*yh + tzn*zh;
    s=s/4.0; 
#ifdef TS_DOUBLE_DOUBLE
    if(ht>=0.0) {
        vtx->curvature=sqrt(h);
    } else {
        vtx->curvature=-sqrt(h);
    }
#endif
#ifdef TS_DOUBLE_FLOAT
    if(ht>=0.0) {
        vtx->curvature=sqrtf(h);
    } else {
        vtx->curvature=-sqrtf(h);
    }
#endif
#ifdef TS_DOUBLE_LONGDOUBLE
    if(ht>=0.0) {
        vtx->curvature=sqrtl(h);
    } else {
        vtx->curvature=-sqrtl(h);
    }
#endif
// c is forced curvature energy for each vertex. Should be set to zero for
// normal circumstances.
/* the following statement is an expression for $\frac{1}{2}\int(c_1+c_2-c_0^\prime)^2\mathrm{d}A$, where $c_0^\prime=2c_0$ (twice the spontaneous curvature)  */
    vtx->energy=0.5*s*(vtx->curvature/s-vtx->c)*(vtx->curvature/s-vtx->c);

    return TS_SUCCESS;
}



ts_bool sweep_attraction_bond_energy(ts_vesicle *vesicle){
	int i;
	for(i=0;i<vesicle->blist->n;i++){
		attraction_bond_energy(vesicle->blist->bond[i], vesicle->tape->w);
	}
	return TS_SUCCESS;
}


inline ts_bool attraction_bond_energy(ts_bond *bond, ts_double w){

	if(fabs(bond->vtx1->c)>1e-16 && fabs(bond->vtx2->c)>1e-16){
		bond->energy=-w;
	}
	else {
		bond->energy=0.0;
	}
	return TS_SUCCESS;
}




ts_double wall_force(ts_vesicle *vesicle, ts_double x1, ts_double y1){
	double r;
	r=sin(vesicle->tape->wall_theta)*(x1-vesicle->tape->wallx)-cos(vesicle->tape->wall_theta)*(y1-vesicle->tape->wally);
	if (r<0){
		return 0;
	}
	else{
		return r;
	}
}


ts_double wall_force_energy(ts_vesicle *vesicle, ts_vertex *vtx, ts_vertex *vtx_old){
	if (vesicle->tape->wall_spring<1e-15){
		return 0;
	}
	double dist, wFx, wFy;
	dist=wall_force(vesicle, vtx->x, vtx->y);
	wFx=-vesicle->tape->wall_spring*dist*sin(vesicle->tape->wall_theta);
	wFy=vesicle->tape->wall_spring*dist*cos(vesicle->tape->wall_theta);

	return -wFx*(vtx->x-vtx_old->x)-wFy*(vtx->y-vtx_old->y);
}


ts_double direct_force_energy(ts_vesicle *vesicle, ts_vertex *vtx, ts_vertex *vtx_old){
	if(fabs(vtx->c)<1e-15) return 0.0;

	if(fabs(vesicle->tape->F)<1e-15) return 0.0;

	ts_double proj, conc, cs=1.0;

	//printf("v--->  %d\n", vtx->idx);

	ts_double norml,ddp=0.0;
	ts_uint i;
	ts_double xnorm=0.0,ynorm=0.0,znorm=0.0;
	/*find normal of the vertex as sum of all the normals of the triangles surrounding it. */
	for(i=0;i<vtx->tristar_no;i++){
		xnorm+=vtx->tristar[i]->xnorm;
		ynorm+=vtx->tristar[i]->ynorm;
		znorm+=vtx->tristar[i]->znorm;
	}
	/*normalize*/
	norml=sqrt(xnorm*xnorm+ynorm*ynorm+znorm*znorm);
	xnorm/=norml;
	ynorm/=norml;
	znorm/=norml;

	// Shubhadeep //
	ts_double F;
	F=vtx->F_vtx;
	
	/*Shubhadeep*/
	// Inhibition force //
	if(vesicle->tape->inhibition_switch){
		//F=sqrt(vtx->Factx*vtx->Factx+vtx->Facty*vtx->Facty+vtx->Factz*vtx->Factz);
		
		ts_double cs, conc;
		cs=vesicle->tape->cs;
		double c0=vesicle->tape->conc0;
		vtx->proj=vesicle->xnorm*(vtx->x-vesicle->cm[0])+vesicle->ynorm*(vtx->y-vesicle->cm[1])+vesicle->znorm*(vtx->z-vesicle->cm[2]);
		conc=c0*vesicle->tape->beta*(vesicle->proj_max-vesicle->proj_min)*vesicle->vmag*exp(-vesicle->tape->beta*vesicle->vmag*(vtx->proj-vesicle->proj_min)/vesicle->tape->D)/(1-exp(-vesicle->tape->beta*vesicle->vmag*(vesicle->proj_max-vesicle->proj_min)/vesicle->tape->D))/vesicle->tape->D;
		conc/=vesicle->volume;
		F=vtx->F_vtx*cs/(conc+cs);
		vtx->cx=conc;
		
	}
	

	vtx->Factx=-F*xnorm;
	vtx->Facty=-F*ynorm;
	vtx->Factz=-F*znorm;
	
	vtx->Fx=-xnorm;
	vtx->Fy=-ynorm;
	vtx->Fz=-znorm;

	/*calculate ddp, perpendicular displacement*/
	ddp=vtx->Factx*(vtx->x-vtx_old->x)+vtx->Facty*(vtx->y-vtx_old->y)+vtx->Factz*(vtx->z-vtx_old->z);
	/*calculate dE*/
	return -ddp;	
	//Shubhadeep //	
	
}

ts_double blow_force_energy(ts_vesicle *vesicle,  ts_vertex *vtx, ts_vertex *vtx_old){
	if(fabs(vesicle->tape->Fblow)<1e-15) return 0.0;
	int i;
	ts_double xnorm=0.0,ynorm=0.0,znorm=0.0, norml;
	/*find normal of the vertex as sum of all the normals of the triangles surrounding it. */
	for(i=0;i<vtx->tristar_no;i++){
		xnorm+=vtx->tristar[i]->xnorm;
		ynorm+=vtx->tristar[i]->ynorm;
		znorm+=vtx->tristar[i]->znorm;
	}
	/*normalize*/
	norml=sqrt(xnorm*xnorm+ynorm*ynorm+znorm*znorm);
	xnorm/=norml;
	ynorm/=norml;
	znorm/=norml;
	double Fp, Ft; 

	double n1, n2, n3, energy; //, Fshy, Fshz;
	n1=0;
	n2=0;
	n3=0;

	
	cross_product(xnorm, ynorm, znorm, 1, 0, 0, &n1, &n2, &n3);
	cross_product(n1, n2, n3, xnorm, ynorm, znorm, &vtx->Fbx, &vtx->Fby, &vtx->Fbz);

	double sigma=100;
	if (xnorm>0 && vtx->x<vesicle->xback+2){
		Fp=vesicle->tape->Fblow*exp(-(vtx->y-vesicle->cm[1])*(vtx->y-vesicle->cm[1])/sigma-(vtx->z-vesicle->tape->z_adhesion)*(vtx->z-vesicle->tape->z_adhesion)/sigma)*xnorm;	
		Ft=vesicle->tape->Fblow*exp(-(vtx->y-vesicle->cm[1])*(vtx->y-vesicle->cm[1])/sigma-(vtx->z-vesicle->tape->z_adhesion)*(vtx->z-vesicle->tape->z_adhesion)/sigma)*sqrt(1-xnorm*xnorm);
	}
	else{
		Fp=0;
		Ft=0;
	}	
	energy=-Ft*(vtx->Fbx*(vtx->x-vtx_old->x)+vtx->Fby*(vtx->y-vtx_old->y)+vtx->Fbz*(vtx->z-vtx_old->z));

	vtx->Fbx=Fp*xnorm;
	vtx->Fby=Fp*ynorm;
	vtx->Fbz=Fp*znorm;

	energy+=-(vtx->Fbx*(vtx->x-vtx_old->x)+vtx->Fby*(vtx->y-vtx_old->y)+vtx->Fbz*(vtx->z-vtx_old->z));
	//double theta=acos(xnorm);
	return  energy;

}


ts_double direct_force_energy_with_Fz_balance(ts_vesicle *vesicle, ts_vertex *vtx, ts_vertex *vtx_old){
	if(fabs(vesicle->tape->F)<1e-15) return 0.0;

	ts_double norml,ddp=0.0;
	ts_uint i;
	ts_double Fz;
	ts_double xnorm=0.0,ynorm=0.0,znorm=0.0;
	force_per_vertex(vesicle, &Fz);
	if(fabs(vtx->c)<1e-15){

	znorm=-Fz;

	/*calculate ddp, perpendicular displacement*/
	ddp=xnorm*(vtx->x-vtx_old->x)+ynorm*(vtx->y-vtx_old->y)+znorm*(vtx->z-vtx_old->z);
	/*calculate dE*/
//	printf("ddp=%e",ddp);
	return vesicle->tape->F*ddp;	

	}

	/*find normal of the vertex as sum of all the normals of the triangles surrounding it. */
	for(i=0;i<vtx->tristar_no;i++){
			xnorm+=vtx->tristar[i]->xnorm;
			ynorm+=vtx->tristar[i]->ynorm;
			znorm+=vtx->tristar[i]->znorm;
	}
	/*normalize*/
	norml=sqrt(xnorm*xnorm+ynorm*ynorm+znorm*znorm);
	xnorm/=norml;
	ynorm/=norml;
	znorm/=norml;
    /*subtract the force per vertex for overall force balance in z-direction*/

	znorm=znorm-Fz;

	/*calculate ddp, perpendicular displacement*/
	ddp=xnorm*(vtx->x-vtx_old->x)+ynorm*(vtx->y-vtx_old->y)+znorm*(vtx->z-vtx_old->z);
	/*calculate dE*/
//	printf("ddp=%e",ddp);
	return vesicle->tape->F*ddp;		
	
}


// Shubhadeep shear force //
//ts_double shear_force_energy(ts_vesicle *vesicle, ts_vertex *vtx, ts_vertex *vtx_old){
//	if(fabs(vesicle->tape->a)<1e-15) return 0.0;
//	return 0.0;
	
/*
	ts_double Fs=vesicle->tape->a*(vtx->z-vesicle->tape->z_adhesion)+vesicle->tape->b;


	


	ts_double nx,ny;
	nx=vesicle->tape->nx;
	ny=vesicle->tape->ny;
	Fs=Fs/sqrt(nx*nx+ny*ny);
	
	vtx->Fshx=Fs*nx;
	vtx->Fshy=Fs*ny;
	vtx->Fshz=0;
	return -vtx->Fshx*(vtx->x-vtx_old->x)-vtx->Fshy*(vtx->y-vtx_old->y);
*/	
//}
/* Shubhadeep Shear force */
void cross_product(double x1, double y1, double z1, double x2, double y2,  double z2, double* xout, double* yout, double* zout){
	*xout=y1*z2-y2*z1;
	*yout=z1*x2-x1*z2;
	*zout=x1*y2-y1*x2;
}


/* Shubhadeep */
/*################################################################ */
ts_double shear_force_energy(ts_vesicle *vesicle, ts_vertex *vtx, ts_vertex *vtx_old){
	//if (vesicle->tape->Fud<1e-15) return 0;
	ts_double xnorm=0;
	ts_double ynorm=0;
	ts_double znorm=0;
	ts_uint i;
	/*find normal of the vertex as sum of all the normals of the triangles surrounding it. */
	for(i=0;i<vtx->tristar_no;i++){
		xnorm+=vtx->tristar[i]->xnorm;
		ynorm+=vtx->tristar[i]->ynorm;
		znorm+=vtx->tristar[i]->znorm;
	}
	/*normalize*/
	ts_double norml;
	norml=sqrt(xnorm*xnorm+ynorm*ynorm+znorm*znorm);
	xnorm/=norml;
	ynorm/=norml;
	znorm/=norml;


	ts_double nx,ny;
	nx=vesicle->tape->nx;
	ny=vesicle->tape->ny;

	nx/=sqrt(nx*nx+ny*ny);
	ny/=sqrt(nx*nx+ny*ny);

	double n1, n2, n3; //, Fshy, Fshz;
	n1=0;
	n2=0;
	n3=0;

	
	
	cross_product(xnorm, ynorm, znorm, nx, ny, 0, &n1, &n2, &n3);
	cross_product(n1, n2, n3, xnorm, ynorm, znorm, &vtx->Fshx, &vtx->Fshy, &vtx->Fshz);
	
	//printf("%1.4f \t %1.4f %1.4f \n",n1, n2, n3);
	ts_double Fs=vesicle->tape->a*(vtx->z-vesicle->tape->z_adhesion)+vesicle->tape->b;
	
	vtx->Fshx=Fs*vtx->Fshx;
	vtx->Fshy=Fs*vtx->Fshy;
	vtx->Fshz=Fs*vtx->Fshz;

	return -vtx->Fshx*(vtx->x-vtx_old->x)-vtx->Fshy*(vtx->y-vtx_old->y)-vtx->Fshz*(vtx->z-vtx_old->z);
	/*
	ts_double projx;
	ts_double projy;
	

	projx=xnorm-(ny*xnorm-nx*ynorm)*ny;
	projy=ynorm+(ny*xnorm-nx*ynorm)*nx;
	


	ts_double proj_norm=sqrt(projx*projx+projy*projy);
	projx/=proj_norm;
	projy/=proj_norm;

	ts_double force=0;

	ts_double Fs=vesicle->tape->a*(vtx->z-vesicle->tape->z_adhesion)+vesicle->tape->b;
	

	if (nx!=0){
		vtx->Fshx=fabs(Fs*nx*xnorm)*nx/fabs(nx);
	}
	if (ny!=0){
		vtx->Fshy=fabs(Fs*ny*ynorm)*ny/fabs(ny);
	}
	*/
	//printf("%1.5e,  %1.5e\n", nx, ny);
	//printf("new %1.5e\n",(fabs(Fs*nx*xnorm)*nx/fabs(nx)));
	//printf("old %1.5e\n",((Fs*nx*xnorm)));
	/*
	if((xnorm*nx+ynorm*ny)>0){

		if (vtx->z-vesicle->tape->z_adhesion<1/5.){
			force=0;
		}
		else{
			force=vesicle->tape->Fud*exp(-(vtx->z-vesicle->tape->z_adhesion));
		}
		vtx->Fshz=force;
		return force*(vtx->z-vtx_old->z)-vtx->Fshx*(vtx->x-vtx_old->x)-vtx->Fshy*(vtx->y-vtx_old->y);
	}
	else{
		if (vtx->z-vesicle->tape->z_adhesion<1/5.){
			force=0;
		}
		else{
			force=vesicle->tape->Fud*exp(-(vtx->z-vesicle->tape->z_adhesion));
        }
        vtx->Fshz=-force;
		return -force*(vtx->z-vtx_old->z)-vtx->Fshx*(vtx->x-vtx_old->x)-vtx->Fshy*(vtx->y-vtx_old->y);
	}
	*/
}
/*################################################################ */

void force_per_vertex(ts_vesicle *vesicle, ts_double *Fz){
	ts_double norml;
	ts_uint i,j;
	ts_double fz=0.0;
	ts_double xnorm,ynorm,znorm;
	/*find normal of the vertex as sum of all the normals of the triangles surrounding it. */
    for (j=0;j<vesicle->vlist->n;j++){
		if(vesicle->vlist->vtx[j]->c > 1e-15){
			xnorm=0.0;ynorm=0.0;znorm=0.0;
			for(i=0;i<vesicle->vlist->vtx[j]->tristar_no;i++){
				xnorm+=vesicle->vlist->vtx[j]->tristar[i]->xnorm;
				ynorm+=vesicle->vlist->vtx[j]->tristar[i]->ynorm;
				znorm+=vesicle->vlist->vtx[j]->tristar[i]->znorm;
			}
			/*normalize*/
			norml=sqrt(xnorm*xnorm+ynorm*ynorm+znorm*znorm);
			znorm/=norml;

			fz+=znorm;
		}
	}
	fz/=vesicle->vlist->n;
	if(fz<0.0000){
		*Fz=fz;
	}
	else if (fz>=0){
		*Fz=0.0;
	}
}

void stretchenergy(ts_vesicle *vesicle, ts_triangle *triangle){
	triangle->energy=vesicle->tape->xkA0/2.0*pow((triangle->area/vesicle->tlist->a0-1.0),2);
}
