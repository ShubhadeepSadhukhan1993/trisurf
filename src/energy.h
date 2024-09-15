/* vim: set ts=4 sts=4 sw=4 noet : */
#ifndef _ENERGY_H
#define _ENERGY_H

ts_bool mean_curvature_and_energy(ts_vesicle *vesicle);
inline ts_bool energy_vertex(ts_vertex *vtx);
inline ts_bool bond_energy(ts_bond *bond,ts_poly *poly);

ts_bool sweep_attraction_bond_energy(ts_vesicle *vesicle);
inline ts_bool attraction_bond_energy(ts_bond *bond, ts_double w);
ts_double direct_force_energy(ts_vesicle *vesicle, ts_vertex *vtx, ts_vertex *vtx_old);
ts_double shear_force_energy(ts_vesicle *vesicle, ts_vertex *vtx, ts_vertex *vtx_old);


ts_double updown_force_energy(ts_vesicle *vesicle, ts_vertex *vtx, ts_vertex *vtx_old);
void stretchenergy(ts_vesicle *vesicle, ts_triangle *triangle);


ts_double direct_force_energy_with_Fz_balance(ts_vesicle *vesicle, ts_vertex *vtx, ts_vertex *vtx_old);
void force_per_vertex(ts_vesicle *vesicle, ts_double *Fz);
void cross_product(double x1, double y1, double z1, double x2, double y2, double z2, double *xout, double *yout, double *zout);

ts_double blow_force_energy(ts_vesicle *vesicle,  ts_vertex *vtx, ts_vertex *vtx_old);
ts_double wall_force_energy(ts_vesicle *vesicle, ts_vertex *vtx, ts_vertex *vtx_old);
ts_double random_force(ts_vesicle *vesicle);
void vesicle_random_force(ts_vesicle *vesicle);
#endif
